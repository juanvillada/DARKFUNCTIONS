# Working with MAGs

We will be working here with all MAGs (≥50% completeness, ≤10% contamination) from co-assemblies #1 and #2. We will start by renaming the contigs to make their names shorter and more informative adding the MAG they belong to.

### Get MAGs and rename contigs

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

# ASSEMBLY 01
for CLUSTER in $(seq 31 100 | awk -v OFS='_' '{print "Cluster", $0}'); do
  for MAG in $(sqlite3 $WORKDIR/04_BINNING/ASSEMBLY_01/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info' | grep MAG); do
    COUNT=$(expr $COUNT + 1)
    NEWMAG=$(printf '%s_%s_%04d\n' heath MAG $COUNT)

    printf '%s\t%s\n' $MAG $NEWMAG >> renamed_MAGs.txt

    anvi-script-reformat-fasta $WORKDIR/04_BINNING/ASSEMBLY_01/CONCOCT_SPLIT/$CLUSTER.SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               --output-file REDUNDANT_MAGS/$NEWMAG.fa \
                               --simplify-names \
                               --prefix $NEWMAG
  done
done

# ASSEMBLY 02
for CLUSTER in $(seq 31 100 | awk -v OFS='_' '{print "Cluster", $0}'); do
  for MAG in $(sqlite3 $WORKDIR/04_BINNING/ASSEMBLY_02/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info' | grep MAG); do
    COUNT=$(expr $COUNT + 1)
    NEWMAG=$(printf '%s_%s_%04d\n' bog MAG $COUNT)

    printf '%s\t%s\n' $MAG $NEWMAG >> renamed_MAGs.txt

    anvi-script-reformat-fasta $WORKDIR/04_BINNING/ASSEMBLY_02/CONCOCT_SPLIT/$CLUSTER.SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               --output-file REDUNDANT_MAGS/$NEWMAG.fa \
                               --simplify-names \
                               --prefix $NEWMAG
  done
done
```

---

Now we will run CHECKM and QUAST to compute some statistics about the MAGs.

### Analyse MAGs with CHECKM

```bash
cd $WORKDIR/05_MAGs

conda activate checkm

checkm lineage_wf REDUNDANT_MAGS \
                  CHECKM \
                  --extension fa \
                  --threads 40 \
                  --pplacer_threads 40

checkm tree_qa CHECKM \
               --file CHECKM/checkm_tree_qa.txt \
               --out_format 1 \
               --tab_table

checkm qa CHECKM/lineage.ms \
          CHECKM \
          --file CHECKM/checkm_qa.txt \
          --out_format 1 \
          --threads 40 \
          --tab_table
```

### Analyse MAGs with QUAST

```bash
cd $WORKDIR/05_MAGs

conda activate quast

quast.py REDUNDANT_MAGS/*.fa \
         -o QUAST \
         -t 40 \
         --no-plots
```

---

Here we will perform a prelimary phylogenomic analysis of the MAGs with GTDB-tk.

### Do phylogenomic analysis with GTDB-tk

```bash
cd $WORKDIR/05_MAGs

conda activate GTDBtk

gtdbtk classify_wf --genome_dir REDUNDANT_MAGS \
                   --out_dir GTDB \
                   --extension fa \
                   --cpus 40
```

---

As MAGs come from two different co-assemblies, it is likely that MAGs corresponding to the same population are present multiple times. To avoid underestimating their abundance and detection across the samples we will identify redundant MAGs and create a final dataset of non-redundant MAGs that will be used for mapping.

### Dereplicate MAGs

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

printf '%s\t%s\n' name path > fastANI_paths.txt

for MAG in $(ls REDUNDANT_MAGS/ | sed 's/.fa//g'); do
  echo $MAG | awk -v OFS='\t' '{print $0, "REDUNDANT_MAGS/" $0 ".fa"}'
done >> fastANI_paths.txt

anvi-dereplicate-genomes --fasta-text-file fastANI_paths.txt \
                         --output-dir FASTANI \
                         --program fastANI \
                         --num-threads 40 \
                         --similarity-threshold 0.90

cat NON_REDUNDANT_MAGS/GENOMES/*.fa > non_redundant_MAGs.fa
```

--

We will now map the quality-filtered metagenomic and metatranscriptomic data to the set of non-redundant MAGs to quantify their abundance across the samples. For metatranscriptomes, we will will work with the data after removing rRNA reads using SORTMERNA.

### Map reads to non-redundant MAGs with BOWTIE

```bash
# Metagenomes
cd $WORKDIR/05_MAGs/MAPPING

module load biokit

# Index contigs
bowtie2-build $WORKDIR/05_MAGs/non_redundant_MAGs.fa \
              non_redundant_MAGs

# Map
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  # Run BOWTIE
  bowtie2 -1 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R1_trimmed.fastq \
          -2 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R2_trimmed.fastq \
          -S $SAMPLE.sam \
          -x non_redundant_MAGs \
          --threads 40 \
          --no-unal

  # Create and index BAM file
  samtools view -F 4 -bS $SAMPLE.sam > $SAMPLE.RAW.bam
  samtools sort $SAMPLE.RAW.bam -o $SAMPLE.bam
  samtools index $SAMPLE.bam

  # Remove SAM and non-indexed BAM files
  rm -f $SAMPLE.sam $SAMPLE.RAW.bam
done

# Metatranscriptomes
cd $WORKDIR/05_MAGs/MAPPING_RNASEQ

module load biokit

# Index contigs
bowtie2-build $WORKDIR/05_MAGs/non_redundant_MAGs.fa \
              non_redundant_MAGs

# Map
for SAMPLE in $(cat $WORKDIR/SAMPLES_RNAseq.txt); do
  # Run BOWTIE
  bowtie2 -U $WORKDIR/01_TRIMMED_DATA_RNASEQ/"$SAMPLE"_trimmed.fasta \
          -S $SAMPLE.sam \
          -x non_redundant_MAGs \
          -f \
          --threads 40 \
          --no-unal

  # Create and index BAM file
  samtools view -F 4 -bS $SAMPLE.sam > $SAMPLE.RAW.bam
  samtools sort $SAMPLE.RAW.bam -o $SAMPLE.bam
  samtools index $SAMPLE.bam

  # Remove SAM and non-indexed BAM files
  rm -f $SAMPLE.sam $SAMPLE.RAW.bam
done
```

---

And now we will import the set of non-redundant MAGs to ANVI'O.

### Import non-redundant MAGs to ANVI'O

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

# Build a contigs database
anvi-gen-contigs-database --contigs-fasta non_redundant_MAGs.fa \
                          --output-db-path CONTIGS.db \
                          --project-name NR_MAGS

# Find single-copy genes with HMMER
anvi-run-hmms --contigs-db CONTIGS.db \
              --num-threads 40

# Get taxonomy for single copy genes
anvi-run-scg-taxonomy --contigs-db CONTIGS.db \
                      --num-parallel-processes 10 \
                      --num-threads 4

# Create map of splits to bins
for SPLIT in $(sqlite3 CONTIGS.db 'select split from splits_basic_info'); do
  MAG=$(echo $SPLIT | awk -F '_' -v OFS='_' '{print $1, $2, $3}')
  printf '%s\t%s\n' $SPLIT $MAG
done > non_redundant_MAGs_splits.txt

# Metagenome profiling

## Build profile databases
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  anvi-profile --input-file MAPPING/$SAMPLE.bam \
               --output-dir PROFILES/$SAMPLE \
               --contigs-db CONTIGS.db \
               --num-threads 40
done

## Merge profiles
anvi-merge PROFILES/*/PROFILE.db \
           --output-dir MERGED_PROFILES \
           --contigs-db CONTIGS.db

## Create collection
anvi-import-collection non_redundant_MAGs_splits.txt \
                       --contigs-db CONTIGS.db \
                       --pan-or-profile-db MERGED_PROFILES/PROFILE.db \
                       --collection-name NR_MAGS

# Metatranscriptome profiling

## Build profile databases
for SAMPLE in $(cat $WORKDIR/SAMPLES_RNAseq.txt); do
  anvi-profile --input-file MAPPING_RNASEQ/$SAMPLE.bam \
               --output-dir PROFILES_RNASEQ/$SAMPLE \
               --contigs-db CONTIGS.db \
               --num-threads 40
done

## Merge profiles
anvi-merge PROFILES_RNASEQ/*/PROFILE.db \
           --output-dir MERGED_PROFILES_RNASEQ \
           --contigs-db CONTIGS.db

## Create collection
anvi-import-collection non_redundant_MAGs_splits.txt \
                       --contigs-db CONTIGS.db \
                       --pan-or-profile-db MERGED_PROFILES_RNASEQ/PROFILE.db \
                       --collection-name NR_MAGS_RNASEQ
```

---

Now we will start annotating the MAGs. We will do this on the gene calls found by PRODIGAL within ANVI'O. We will annotate them against the KEGG database using DIAMOND and against the SEED database using myRAST.

### Get sequences for gene calls

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

anvi-export-gene-calls --contigs-db CONTIGS.db \
                       --output-file gene_calls.txt \
                       --gene-caller prodigal

anvi-get-sequences-for-gene-calls --contigs-db CONTIGS.db \
                                  --output-file gene_calls.fa

anvi-get-sequences-for-gene-calls --contigs-db CONTIGS.db \
                                  --output-file gene_calls.faa \
                                  --get-aa-sequences
```

### Annotate genes against COG with BLAST

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

anvi-run-ncbi-cogs --contigs-db CONTIGS.db \
                   --num-threads 40
```

### Annotate genes against KEGG with DIAMOND

```bash
cd $WORKDIR/05_MAGs

conda activate KEGGtools

# Run DIAMOND
diamond blastx --query gene_calls.fa \
               --out gene_calls_KEGG.txt \
               --db $KEGG/PROKARYOTES \
               --outfmt 6 \
               --evalue 0.00001 \
               --id 60 \
               --max-target-seqs 1 \
               --max-hsps 1 \
               --threads 40

# Run KEGG-tools
KEGG-tools-assign.py -i gene_calls_KEGG.txt \
                     -a $KEGG

# Create an annotation file and import to ANVI'O
conda activate anvio-master

printf '%s\t%s\t%s\t%s\t%s\n' gene_callers_id source accession function e_value > gene_calls_KEGG_annot.txt

sed '1d' gene_calls_KEGG_KOtable.txt |
awk -v OFS='\t' -F '\t' '{print $1, "KEGG", $3, $6, $5}'  >> gene_calls_KEGG_annot.txt

anvi-import-functions --contigs-db CONTIGS.db \
                      --input-files gene_calls_KEGG_annot.txt
```

### Map genes to KEGG modules

```bash
cd $WORKDIR/05_MAGs

conda activate KEGGtools

# Split KEGG annotation by MAG
cut -f 1,2 gene_calls.txt |
cut -f 1-3 -d '_' |
sed '1d' |
sort -k 1,1 > KEGG-MODULES/gene_calls_KEGG_by_MAG.txt

for MAG in $(cat non_redundant_MAGs.txt); do
  grep $MAG KEGG-MODULES/gene_calls_KEGG_by_MAG.txt |
  cut -f 1 |
  join -1 1 -2 1 -t $'\t' - <(sort -k 1,1 gene_calls_KEGG.txt) > KEGG-MODULES/$MAG.txt

  KEGG-tools-assign.py -i KEGG-MODULES/$MAG.txt \
                       -o KEGG-MODULES \
                       -a $KEGG \
                       --summarise
done
```

### Summarize refined MAGs

```bash
cd $WORKDIR/05_MAGs

conda activate anvio-master

anvi-summarize --contigs-db CONTIGS.db \
               --output-dir REFINED_MAGS_SUMMARY \
               --pan-or-profile-db MERGED_PROFILES/PROFILE.db \
               --collection-name NR_MAGS \
               --init-gene-coverages

anvi-summarize --contigs-db CONTIGS.db \
               --output-dir REFINED_MAGS_RNASEQ_SUMMARY \
               --pan-or-profile-db MERGED_PROFILES_RNASEQ/PROFILE.db \
               --collection-name NR_MAGS_RNASEQ \
               --init-gene-coverages
```
