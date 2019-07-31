# Read-based analyses

First we will use METAXA to look for 16S rRNA gene fragments and obtain a taxonomic profile for each sample. Next we will run DIAMOND to map the reads against the KEGG database and obtain a functional profile for each sample. We will then use an in-house script (https://github.com/igorspp/KEGG-tools) to parse the DIAMOND results. It (i) assigns a KO identifier to each hit, (ii) runs MinPath to remove spurious pathways, and (iii) summarises the abundance of each pathway. We will work here with the dataset which has been resampled to an even number of reads.

### Annotate reads with METAXA

```bash
cd $WORKDIR/02_METAXA

module load biokit

# Run METAXA
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  metaxa2 -1 $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R1_trimmed.fastq \
          -2 $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R2_trimmed.fastq \
          -o $SAMPLE \
          --align none \
          --graphical F \
          --cpu 4 \
          --plus

  metaxa2_ttt -i "$SAMPLE".taxonomy.txt \
              -o $SAMPLE
done

# Summarise results
for LEVEL in 1 2 3 4 5 6 7 8; do
  metaxa2_dc *level_"$LEVEL".txt \
             -o summary_level_"$LEVEL".txt
done
```

### Annotate reads against the KEGG database with DIAMOND

```bash
cd $WORKDIR/02_KEGG

module load biokit/4.9.3
module load biopython-env/2.7.13

# Convert to FASTA
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  seqtk seq -A $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R1_trimmed.fastq | awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R1", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta
  seqtk seq -A $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R2_trimmed.fastq | awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R2", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta
done

# Run DIAMOND
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  diamond blastx --query "$SAMPLE"_trimmed.fasta \
                 --out "$SAMPLE".txt \
                 --db $KEGG/PROKARYOTES \
                 --outfmt 6 \
                 --evalue 0.00001 \
                 --id 60 \
                 --max-target-seqs 1 \
                 --max-hsps 1 \
                 --threads 4
done

# Run KEGG-tools
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  KEGG-tools-assign.py -i "$SAMPLE".txt \
                       -p "$SAMPLE" \
                       -a $KEGG \
                       --summarise \
                       --minpath
done

KEGG-tools-collect.py -s $WORKDIR/SAMPLES.txt
```
