# Amplicon sequencing analyses

```bash
# Quality check of raw reads
conda activate fastqc

fastqc 01_RAW_DATA/*.fastq.gz \
       -o 01_RAW_DATA/FASTQC \
       -t 4

multiqc 01_RAW_DATA/FASTQC \
        -o 01_RAW_DATA/MULTIQC \
        --interactive

# Merge paired-end reads
for SAMPLE in $(cat SAMPLES.txt); do
  gunzip -c 01_RAW_DATA/"$SAMPLE"_*_R1_*.fastq.gz > "$SAMPLE"_R1.fastq
  gunzip -c 01_RAW_DATA/"$SAMPLE"_*_R2_*.fastq.gz > "$SAMPLE"_R2.fastq

  usearch -fastq_mergepairs "$SAMPLE"_R1.fastq \
          -reverse "$SAMPLE"_R2.fastq \
          -fastqout 02_MERGED_DATA/$SAMPLE.fastq \
          -fastq_maxdiffs 10 \
          -threads 4 &> 02_MERGED_DATA/"$SAMPLE".log

  rm -f "$SAMPLE"_R1.fastq "$SAMPLE"_R2.fastq
done

# Trim primers
conda activate cutadapt

## ARCHAEAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "A-"); do
  cutadapt 02_MERGED_DATA/$SAMPLE.fastq \
           -o 03_TRIMMED_DATA/$SAMPLE.fastq \
           -a ^GNGCANCAGNCGNGAAN...CAGCNGCCGCGGTAA$ \
           -j 4 \
           --discard-untrimmed > 03_TRIMMED_DATA/$SAMPLE.log
done

## BACTERIAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "B-"); do
  cutadapt 02_MERGED_DATA/$SAMPLE.fastq \
           -o 03_TRIMMED_DATA/$SAMPLE.fastq \
           -a ^GTGCCAGCMGCCGCGGTAA...ATTAGANACCCNNGTAGTCC$ \
           -j 4 \
           --discard-untrimmed > 03_TRIMMED_DATA/$SAMPLE.log
done

# Quality check of trimmed reads
conda activate fastqc

fastqc 03_TRIMMED_DATA/*.fastq \
       -o 03_TRIMMED_DATA/FASTQC \
       -t 4

multiqc 03_TRIMMED_DATA/FASTQC \
        -o 03_TRIMMED_DATA/MULTIQC \
        --interactive

# Quality filtering
for SAMPLE in $(cat SAMPLES.txt); do
  usearch -fastq_filter 03_TRIMMED_DATA/$SAMPLE.fastq \
          -fastaout 04_FILTERED_DATA/$SAMPLE.fasta \
          -fastq_maxee 1 \
          -threads 4 &> 04_FILTERED_DATA/$SAMPLE.log
done

# Pool samples
for SAMPLE in $(cat SAMPLES.txt); do
  cat 04_FILTERED_DATA/$SAMPLE.fasta |
  awk -v NAME=$SAMPLE '/^>/{print ">" NAME "_read" ++i ";barcodelabel=" NAME; next}{print}'
done > filtered.fasta

# Dereplicate reads
usearch -fastx_uniques filtered.fasta \
        -fastaout uniques.fasta \
        -relabel Uniq \
        -sizeout \
        -threads 4

# Cluster OTUs/ASVs
usearch -cluster_otus uniques.fasta \
        -otus 05_OTUs/OTUs.fasta \
        -relabel OTU

usearch -unoise3 uniques.fasta \
        -zotus 05_OTUs/ASVs.fasta

# Assign taxonomy with RDP in QIIME
conda activate qiime-1.9.1

## GREENGENES
for DATA in OTUs ASVs; do
  assign_taxonomy.py -i 05_OTUs/$DATA.fasta \
                     -o 05_OTUs/GREENGENES \
                     -m rdp
done

## SILVA
for DATA in OTUs ASVs; do
  assign_taxonomy.py -i 05_OTUs/$DATA.fasta \
                     -o 05_OTUs/SILVA \
                     -t 00_SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt \
                     -r 00_SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
                     -m rdp \
                     --rdp_max_memory 128000
done

# Remove mitochondria and chloroplast
for DATA in OTUs ASVs; do
  cat 05_OTUs/GREENGENES/"$DATA"_tax_assignments.txt |
  grep -v "mitochondria" |
  grep -v "Chloroplast" |
  cut -f 1 > 05_OTUs/$DATA.good.txt

  usearch -fastx_getseqs 05_OTUs/$DATA.fasta \
          -labels 05_OTUs/$DATA.good.txt \
          -fastaout 05_OTUs/$DATA.good.fasta
done

# Make OTU/ASV tables
usearch -otutab filtered.fasta \
        -otus 05_OTUs/OTUs.good.fasta \
        -otutabout 05_OTUs/OTU_table.txt \
        -threads 2

usearch -otutab filtered.fasta \
        -zotus 05_OTUs/ASVs.good.fasta \
        -otutabout 05_OTUs/_ASV_table.txt \
        -threads 2
done
