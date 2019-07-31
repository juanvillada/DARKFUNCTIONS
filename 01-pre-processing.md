# Pre-processing of raw data

We will start with some basic quality checking of the sequencing data using FASTQC. Next we will remove the sequencing adaptors and perform some quality filtering using CUTADAPT, and run FASTQC again to check the quality of the trimmed data. This will be done for each sequencing run separately by assigning the variable $RUN, which we will be merged afterwards. Finally, since there is a big discrepancy in the number of reads obtained for each sample, we will resample each dataset to the same number of reads (2,000,000) for the read-based analyses.

### Set working directory

```bash
export WORKDIR=/wrk/stelmach/DONOTREMOVE/DARKFUNCTIONS
```

### Define sequencing run

```bash
RUN=1STRUN #OR
RUN=2NDRUN #OR
RUN=3RDRUN
```

### Create list of sample names

```bash
cd $WORKDIR

ls 01_RAW_DATA/$RUN/*.fastq.gz |
egrep -o '[mo][0-9]+' |
sort |
uniq > SAMPLES_"$RUN".txt
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

for SAMPLE in $(cat $WORKDIR/SAMPLES_"$RUN".txt); do
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

### Pool data from all sequencing runs

```bash
# Create list of sample names
cd $WORKDIR/

cat samples_1STRUN.txt samples_2NDRUN.txt samples_3RDRUN.txt |
sort |
uniq > SAMPLES.txt

# Pool reads
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
