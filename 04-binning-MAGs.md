# Binning of MAGs

## Co-assemblies

### Define assembly and create list of sample names

```bash
ASSEMBLY=01
SAMPLES=$(cat $WORKDIR/SAMPLES.txt | egrep -v [mo]122)

#OR

ASSEMBLY=02
SAMPLES=$(cat $WORKDIR/SAMPLES.txt | egrep [mo]122)
```

### Bin MAGs with ANVI’O

Here we will start binning contigs into MAGs with ANVI’O. First we rename the contigs (which have long names given by MEGAHIT) and select only those >2,500 bp. Then we create an ANVI’O contigs database which computes k-mer frequencies for each contig, soft-splits contigs >20,000 bp and identify ORFs using Prodigal. We then run HMMER to identify rRNAs and single-copy core genes. ORFs found by Prodigal are exported and annotated taxonomically using CENTRIFUGE. Reads are then mapped to the contigs using BOWTIE, and profile databases are created which contains sample-specific information about the contigs.
Profiles from each sample are merged into one single file and CONCOCT is used to pre-cluster a fixed number of clusters to reduce the size of the files and allow binning.

```bash
cd $WORKDIR/04_BINNING/BINNING_"$ASSEMBLY"

conda activate anvio6

# Rename contigs and select those >2,500 bp
anvi-script-reformat-fasta $WORKDIR/03_ASSEMBLIES/ASSEMBLY_"$ASSEMBLY"/final.contigs.fa \
                           -o CONTIGS_2500nt.fa \
                           -r CONTIGS_reformat.txt \
                           -l 2500 \
                           --prefix assembly_"$ASSEMBLY" \
                           --simplify-names

# Build a contigs database
anvi-gen-contigs-database -f CONTIGS_2500nt.fa \
                          -o CONTIGS.db \
                          -n assembly_"$ASSEMBLY"

# Find single-copy genes with HMMER
anvi-run-hmms -c CONTIGS.db \
              -T 40

# Get sequences for gene calls
anvi-get-sequences-for-gene-calls -c CONTIGS.db \
                                  -o gene_calls.fa

# Annotate gene calls with CENTRIFUGE
centrifuge -f gene_calls.fa \
           -S centrifuge_hits.tsv \
           -x $CENTRIFUGE_BASE/p+h+v \
           -p 16

# Import CENTRIFUGE results
anvi-import-taxonomy-for-genes -i centrifuge_report.tsv centrifuge_hits.tsv \
                               -c CONTIGS.db \
                               -p centrifuge

# Map reads with BOWTIE

## Index contigs
bowtie2-build CONTIGS_2500nt.fa \
              MAPPING/contigs

## Map
for SAMPLE in $SAMPLES; do
  # Run BOWTIE
  bowtie2 -1 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R1_trimmed.fastq \
          -2 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R2_trimmed.fastq  \
          -S MAPPING/$SAMPLE.sam \
          -x MAPPING/contigs \
          --threads 8 \
          --no-unal

  # Create and index BAM file
  samtools view -F 4 -bS MAPPING/$SAMPLE.sam > MAPPING/$SAMPLE.RAW.bam
  samtools sort MAPPING/$SAMPLE.RAW.bam -o MAPPING/$SAMPLE.bam
  samtools index MAPPING/$SAMPLE.bam

  # Remove SAM and unsorted BAM files
  rm -f MAPPING/$SAMPLE.sam MAPPING/$SAMPLE.RAW.bam
done

# Build profile databases
for SAMPLE in $SAMPLES; do
  anvi-profile -i MAPPING/$SAMPLE.bam \
               -o PROFILES/$SAMPLE \
               -c CONTIGS.db \
               -T 24 \
               --skip-SNV-profiling
done

# Merge profiles
anvi-merge PROFILES/*/PROFILE.db \
           -o MERGED_PROFILES \
           -c CONTIGS.db \
           --skip-hierarchical-clustering \
           --skip-concoct-binning

# Get taxonomy for single copy genes
anvi-run-scg-taxonomy -c CONTIGS.db \
                      -P 10 \
                      -T 4

# Pre-cluster with CONCOCT
anvi-cluster-contigs -c CONTIGS.db \
                     -p MERGED_PROFILES/PROFILE.db \
                     -C CONCOCT \
                     -T 4 \
                     --driver concoct \
                     --clusters 100

# Summarize CONCOCT bins
anvi-summarize -c CONTIGS.db \
               -p MERGED_PROFILES/PROFILE.db \
               -C CONCOCT \
               -o CONCOCT_SUMMARY \
               --quick-summary

# Create individual contigs and profile databases
anvi-split -p PROFILE.db \
           -c CONTIGS.db \
           -C CONCOCT \
           -o CONCOCT_SPLIT \
           --skip-variability-tables

# Bin MAGs
for CLUSTER in $(seq 51 100 | awk -v OFS='_' '{print "Bin", $0}'); do
  anvi-refine -p CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
              -c CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
              -b ALL_SPLITS \
              -C DEFAULT \
              --server-only

  anvi-summarize -p CONCOCT_SPLIT/$CLUSTER/PROFILE.db \
                 -c CONCOCT_SPLIT/$CLUSTER/CONTIGS.db \
                 -o CONCOCT_SPLIT/$CLUSTER.SUMMARY \
                 -C DEFAULT &> /dev/null
done
```
