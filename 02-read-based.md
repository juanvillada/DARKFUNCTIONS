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

# Run DIAMOND
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  diamond blastx --query $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_trimmed.fasta \
                 --out "$SAMPLE".txt \
                 --db $KEGG/KEGG \
                 --outfmt 6 \
                 --evalue 0.00001 \
                 --id 60 \
                 --max-target-seqs 1 \
                 --max-hsps 1 \
                 --threads 4
done

# Run KEGG-tools
parallel -j 4 -a $WORKDIR/SAMPLES.txt --gnu "KEGG-tools assign -i {}.txt -p {} -k $KEGG;
                                             KEGG-tools expand -i {}_KOtable.txt -p {} -k $KEGG"

KEGG-tools summarise -s $WORKDIR/SAMPLES.txt \
                     -p summary
```
