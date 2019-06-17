# Pre-processing of raw data

We will start with some basic quality checking of the sequencing data using FASTQC. Next we will remove the sequencing adaptors and perform some quality filtering using CUTADAPT, and run FASTQC again to check the quality of the trimmed data. This will be done for each sequencing run separately.

### Set working directory

```bash
export WORKDIR=/wrk/stelmach/DONOTREMOVE/DARKFUNCTIONS
```

### Create list of sample names

```bash
cd $WORKDIR

ls 01_RAW_DATA/$RUN/*.fastq.gz |
egrep -o [mo][0-9]+ |
sort |
uniq > SAMPLES_"$RUN".txt
done
```

### Check raw data with FASTQC

```bash
cd $WORKDIR/01_RAW_DATA

module load bioconda/3
source activate QC_env

fastqc $RUN/*.fastq.gz \
       -o $RUN/FASTQC \
       -t 4

multiqc $RUN/FASTQC \
        -o $RUN/MULTIQC \
        --interactive
```

### Trim adaptors and do quality filtering with CUTADAPT

```bash
cd $WORKDIR/01_TRIMMED_DATA

module load bioconda/3
source activate cutadapt

for SAMPLE in $(cat $WORKDIR/samples_"$RUN".txt); do
  cutadapt $WORKDIR/01_RAW_DATA/"$RUN"/*"$SAMPLE"_*R1*.fastq.gz \
           $WORKDIR/01_RAW_DATA/"$RUN"/*"$SAMPLE"_*R2*.fastq.gz \
           -o "$RUN"/"$SAMPLE"_R1_trimmed.fastq \
           -p "$RUN"/"$SAMPLE"_R2_trimmed.fastq \
           -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
           -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
           -q 20 \
           -m 50 \
           -j 4 > "$RUN"/"$SAMPLE"_trim.log
done
```

### Check trimmed data with FASTQC

```bash
cd $WORKDIR/01_TRIMMED_DATA

module load bioconda/3
source activate QC_env

fastqc $RUN/*_trimmed.fastq \
       -o $RUN/FASTQC \
       -t 4

multiqc $RUN/FASTQC \
        -o $RUN/MULTIQC \
        --interactive
done
```

---

Now we will merge the trimmed reads from the three sequencing runs. And for the read-based analyses, since there is a big discrepancy in the number of reads obtained for each sample, we will resample each dataset to the same number of reads (2,000,000) and then convert the FASTQ files to FASTA for use with DIAMOND.

### Pool sequencing runs

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

### Resample FASTQ files and convert to FASTA

```bash
cd $WORKDIR/01_TRIMMED_DATA_SUB

module load biokit

for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  seqtk sample -s100 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R1_trimmed.fastq 2000000 > "$SAMPLE"_R1_trimmed.fastq
  seqtk sample -s100 $WORKDIR/01_TRIMMED_DATA/"$SAMPLE"_R2_trimmed.fastq 2000000 > "$SAMPLE"_R2_trimmed.fastq
done

for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  fastq_to_fasta_fast "$SAMPLE"_R1_trimmed.fastq | awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R1", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta
  fastq_to_fasta_fast "$SAMPLE"_R2_trimmed.fastq | awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R2", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta
done
```
