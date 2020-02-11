# Working with MAGs

We will be working here with all MAGs (≥50% completeness, ≤10% contamination) from co-assemblies #1 and #2. We will start by renaming the contigs to make their names shorter and more informative adding the MAG they belong to.

### Get MAGs and rename contigs

```bash
cd $WORKDIR/05_MAGs

conda activate anvio6

# ASSEMBLY 01
for CLUSTER in $(seq 31 100 | awk -v OFS='_' '{print "Cluster", $0}'); do
  for MAG in $(sqlite3 $WORKDIR/04_BINNING/ASSEMBLY_01/CONCOCT_SPLIT/$CLUSTER/PROFILE.db 'SELECT bin_name FROM collections_bins_info' | grep MAG); do
    COUNT=$(expr $COUNT + 1)
    NEWMAG=$(printf '%s_%s_%04d\n' heath MAG $COUNT)

    printf '%s\t%s\n' $MAG $NEWMAG >> renamed_MAGs.txt

    anvi-script-reformat-fasta $WORKDIR/04_BINNING/ASSEMBLY_01/CONCOCT_SPLIT/$CLUSTER.SUMMARY/bin_by_bin/$MAG/$MAG-contigs.fa \
                               -o REDUNDANT_MAGS/$NEWMAG.fa \
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
anvi-activate-master

printf '%s\t%s\n' name path > fastANI_paths.txt

for MAG in $(ls REDUNDANT_MAGS/ | sed 's/.fa//g'); do
  echo $MAG | awk -v OFS='\t' '{print $0, "REDUNDANT_MAGS/" $0 ".fa"}'
done >> fasta_paths.txt

anvi-dereplicate-genomes --fasta-text-file fasta_paths.txt \
                         --output-dir FASTANI \
                         --program fastANI \
                         --num-threads 40 \
                         --similarity-threshold 0.90

mv FASTANI/GENOMES NON_REDUNDANT_MAGS

cat NON_REDUNDANT_MAGS/*.fa > non_redundant_MAGs.fa
```

--

We will now map the quality-filtered metagenomic and metatranscriptomic data to the set of non-redundant MAGs to quantify their abundance across the samples.

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
