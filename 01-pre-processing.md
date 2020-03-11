# Pre-processing of raw data

We will start with some basic quality checking of the sequencing data using FASTQC. Next we will remove the sequencing adaptors and perform some quality filtering using CUTADAPT, and run FASTQC again to check the quality of the trimmed data. This will be done for each sequencing run separately by assigning the variable $RUN, which we will be merged afterwards. For the read-based analyses, we will resample each dataset to the same number of reads (2,000,000) since there is a big discrepancy in the number of reads obtained for each sample. Finally, we will run SORTMERNA to remove rRNA sequences from the metatranscriptomes, for later mapping against the MAGs.

### Set working directory

```bash
export WORKDIR=/scratch/project_2000577
```

### Define sequencing run and create list of sample names

```bash
cd $WORKDIR

RUN=1STRUN #OR
RUN=2NDRUN #OR
RUN=3RDRUN

SAMPLES=$(ls 01_RAW_DATA/$RUN/*.fastq.gz | egrep -o '[mo][0-9]+' | sort | uniq)
```

### Check raw data with FASTQC

```bash
cd $WORKDIR/01_RAW_DATA

module load bioconda/3
source activate QC_env

fastqc $RUN/*.fastq.gz \
       -o $RUN/FASTQC \
       -t 8

multiqc $RUN/FASTQC \
        -o $RUN/MULTIQC \
        --interactive
```

### Trim adaptors and do quality filtering with CUTADAPT

```bash
cd $WORKDIR/01_TRIMMED_DATA

module load bioconda/3
source activate cutadapt

for SAMPLE in $SAMPLES; do
  cutadapt $WORKDIR/01_RAW_DATA/$RUN/*"$SAMPLE"_*R1*.fastq.gz \
           $WORKDIR/01_RAW_DATA/$RUN/*"$SAMPLE"_*R2*.fastq.gz \
           -o $RUN/"$SAMPLE"_R1_trimmed.fastq \
           -p $RUN/"$SAMPLE"_R2_trimmed.fastq \
           -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
           -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
           -m 50 \
           -j 4 \
           --nextseq-trim=20 > $RUN/"$SAMPLE"_trim.log
done
```

### Check trimmed data with FASTQC

```bash
cd $WORKDIR/01_TRIMMED_DATA

module load bioconda/3
source activate QC_env

fastqc $RUN/*_trimmed.fastq \
       -o $RUN/FASTQC \
       -t 8

multiqc $RUN/FASTQC \
        -o $RUN/MULTIQC \
        --interactive
```

### Create list of sample names

```bash
cd $WORKDIR

ls 01_TRIMMED_DATA/*/*.fastq | egrep -o '[mo][0-9]+' | sort | uniq > SAMPLES.txt
```

### Pool data from all sequencing runs

```bash
cd $WORKDIR/01_TRIMMED_DATA

for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  cat */"$SAMPLE"_R1_trimmed.fastq > "$SAMPLE"_R1_trimmed.fastq
  cat */"$SAMPLE"_R2_trimmed.fastq > "$SAMPLE"_R2_trimmed.fastq
done
```

### Resample FASTQ files

```bash
cd $WORKDIR/01_TRIMMED_DATA_SUB

module load biokit

for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  seqtk sample -s100 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R1_trimmed.fastq 2000000 > "$SAMPLE"_R1_trimmed.fastq
  seqtk sample -s100 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R2_trimmed.fastq 2000000 > "$SAMPLE"_R2_trimmed.fastq
done
```

### Remove rRNA reads using SORTMERNA

```bash
cd $WORKDIR/01_TRIMMED_DATA_RNASEQ

SORTMEDB="$HOME/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,$HOME/sortmerna-2.1b/index/silva-bac-16s-db:$HOME/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,$HOME/sortmerna-2.1b/index/silva-arc-16s-db:$HOME/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,$HOME/sortmerna-2.1b/index/silva-euk-18s-db"

# Index the rRNA databases
indexdb_rna -v --ref $SORTMEDB

# Find rRNA sequences
for SAMPLE in $(cat $WORKDIR/SAMPLES_RNAseq.txt); do
  sortmerna --reads "$SAMPLE"_trimmed.fasta
            --ref $SORTMEDB \
            --aligned "$SAMPLE"_rRNA \
            --other "$SAMPLE"_norRNA \
            --fastx \
            -a 20 \
            -v
done
```
