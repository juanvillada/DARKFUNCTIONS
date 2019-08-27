# Amplicon sequencing analyses

```bash
# Quality check of raw reads
module load bioconda/3
source activate QC_env

fastqc 01_RAW/*.fastq.gz \
       -o 01_RAW/FASTQC \
       -t 4

multiqc 01_RAW/FASTQC \
        -o 01_RAW/MULTIQC \
        --interactive

# Merge paired-end reads
for SAMPLE in $(cat SAMPLES.txt); do
  gunzip -c 01_RAW/"$SAMPLE"_*_R1_*.fastq.gz > "$SAMPLE"_R1.fastq
  gunzip -c 01_RAW/"$SAMPLE"_*_R2_*.fastq.gz > "$SAMPLE"_R2.fastq

  usearch -fastq_mergepairs "$SAMPLE"_R1.fastq \
          -reverse "$SAMPLE"_R2.fastq \
          -fastqout 02_MERGED/"$SAMPLE".fastq \
          -fastq_maxdiffs 10 \
          -threads 2 &> 02_MERGED/"$SAMPLE".log

  rm -f "$SAMPLE"_R1.fastq "$SAMPLE"_R2.fastq
done

# Trim primers
module load bioconda/3
source activate cutadapt

## ARCHAEAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "A-"); do
  cutadapt 02_MERGED/"$SAMPLE".fastq \
           -o 03_TRIMMED/"$SAMPLE".fastq \
           -a ^GNGCANCAGNCGNGAAN...CAGCNGCCGCGGTAA$ \
           -j 4 \
           --discard-untrimmed > 03_TRIMMED/"$SAMPLE".log
done

## BACTERIAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "B-"); do
  cutadapt 02_MERGED/"$SAMPLE".fastq \
           -o 03_TRIMMED/"$SAMPLE".fastq \
           -a ^GTGCCAGCMGCCGCGGTAA...ATTAGANACCCNNGTAGTCC$ \
           -j 4 \
           --discard-untrimmed > 03_TRIMMED/"$SAMPLE".log
done

# Quality check of trimmed reads
module load bioconda/3
source activate QC_env

fastqc 03_TRIMMED/*.fastq \
       -o 03_TRIMMED/FASTQC \
       -t 4

multiqc 03_TRIMMED/FASTQC \
        -o 03_TRIMMED/MULTIQC \
        --interactive

# Quality filtering
for SAMPLE in $(cat SAMPLES.txt); do
  usearch -fastq_filter 03_TRIMMED/"$SAMPLE".fastq \
          -fastaout 04_FILTERED/"$SAMPLE".fasta \
          -fastq_maxee 1 \
          -threads 2 &> 04_FILTERED/"$SAMPLE".log
done

# Pool samples

## ARCHAEAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "A-"); do
  cat 04_FILTERED/"$SAMPLE".fasta |
  awk -v NAME=$SAMPLE '/^>/{print ">" NAME "_read" ++i ";barcodelabel=" NAME; next}{print}' >> 04_FILTERED/ARCHAEA_filtered.fasta
done

## BACTERIAL PRIMERS
for SAMPLE in $(cat SAMPLES.txt | grep "B-"); do
  cat 04_FILTERED/"$SAMPLE".fasta |
  awk -v NAME=$SAMPLE '/^>/{print ">" NAME "_read" ++i ";barcodelabel=" NAME; next}{print}' >> 04_FILTERED/BACTERIA_filtered.fasta
done

# Dereplicate reads
for DATA in ARCHAEA BACTERIA; do
  usearch -fastx_uniques 04_FILTERED/"$DATA"_filtered.fasta \
          -fastaout 05_OTUs/"$DATA"_uniques.fasta \
          -relabel Uniq \
          -sizeout \
          -threads 1
done

# Cluster OTUs/ASVs
for DATA in ARCHAEA BACTERIA; do
  usearch -cluster_otus 05_OTUs/"$DATA"_uniques.fasta \
          -otus 05_OTUs/"$DATA"_OTUs.fasta \
          -relabel "$DATA"_OTU

  usearch -unoise3 05_OTUs/"$DATA"_uniques.fasta \
          -zotus 05_OTUs/"$DATA"_ASVs.fasta

  sed -i 's/Zotu/ASV/g' 05_OTUs/"$DATA"_ASVs.fasta
done

# Assign taxonomy with RDP in QIIME
module load qiime/1.9.1

## GREENGENES
for DATA in ARCHAEA BACTERIA; do
  for TYPE in OTUs ASVs; do
    assign_taxonomy.py -i 05_OTUs/"$DATA"_"$TYPE".fasta \
                       -m rdp \
                       -o GREENGENES
  done
done

## SILVA
for DATA in ARCHAEA BACTERIA; do
  for TYPE in OTUs ASVs; do
    assign_taxonomy.py -i 05_OTUs/"$DATA"_"$TYPE".fasta \
                       -t $WRKDIR/DONOTREMOVE/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt \
                       -r $WRKDIR/DONOTREMOVE/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
                       -m rdp \
                       -o SILVA \
                       --rdp_max_memory 128000
  done
done

# Remove mitochondria and chloroplast
for DATA in ARCHAEA BACTERIA; do
  for TYPE in OTUs ASVs; do
    usearch -fastx_getseqs 05_OTUs/"$DATA"_"$TYPE".fasta \
            -labels 05_OTUs/"$DATA"_"$TYPE"_final.txt \
            -fastaout 05_OTUs/"$DATA"_"$TYPE"_final.fasta
  done
done

# Make OTU/ASV tables
for DATA in ARCHAEA BACTERIA; do
  usearch -otutab 04_FILTERED/"$DATA"_filtered.fasta \
          -otus 05_OTUs/"$DATA"_OTUs_final.fasta \
          -otutabout 06_OTU_tables/"$DATA"_OTU_table.txt \
          -threads 2

  usearch -otutab 04_FILTERED/"$DATA"_filtered.fasta \
          -zotus 05_OTUs/"$DATA"_ASVs_final.fasta \
          -otutabout 06_OTU_tables/"$DATA"_ASV_table.txt \
          -threads 2
done
