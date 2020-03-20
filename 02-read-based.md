# Read-based analyses

First we will use METAXA to look for 16S rRNA gene sequences and obtain a taxonomic profile for each sample. Next we will run DIAMOND to map the reads against the KEGG database and generate functional profiles. We will then use an in-house script (https://github.com/igorspp/KEGG-tools) to parse the DIAMOND results. It (i) assigns a KO identifier to each hit, (ii) runs MinPath to remove spurious pathways, and (iii) summarises the abundance of each pathway. We will work here with the dataset which has been resampled to 2,000,000 reads.

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

module load biokit

# Convert reads to FASTA and rename headers
for SAMPLE in $(cat $WORKDIR/SAMPLES.txt); do
  seqtk seq -A $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R1_trimmed.fastq |
  awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R1", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta

  seqtk seq -A $WORKDIR/01_TRIMMED_DATA_SUB/"$SAMPLE"_R2_trimmed.fastq |
  awk -v SAMPLE=$SAMPLE -v OFS='-' '/^>/{print ">" SAMPLE, "R2", "READ", ++i; next}{print}' >> "$SAMPLE"_trimmed.fasta
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
```

```bash
conda activate keggR

# Parse results in R with keggR
library("tidyverse")
library("keggR")

## Create list of samples
SAMPLES <- read_lines("../SAMPLES.txt")

## Load KEGG auxiliary files
loadKEGG("/projappl/project_2000577/KEGG")

## Read BLAST results
data <- lapply(SAMPLES, function(SAMPLE) {
  readBlast(paste(SAMPLE, ".txt", sep = ""))
}) %>% set_names(SAMPLES)

## Assign KO and summarise modules
data <- lapply(data, function(SAMPLE) {
  SAMPLE %>%
    assignKEGG %>%
    runMinpath %>%
    summariseKEGG
})

## Merge summaries
summary <- data %>%
  mergeSummaries

## Write results
save.image("keggR.RData")

summary %>%
  getSummary("modules", "level1") %>%
  write_delim("summaries_modules_level1.txt", delim = "\t")

summary %>%
  getSummary("modules", "level2") %>%
  write_delim("summaries_modules_level2.txt", delim = "\t")

summary %>%
  getSummary("modules", "level3") %>%
  write_delim("summaries_modules_level3.txt", delim = "\t")

summary %>%
  getSummary("modules", "level4") %>%
  write_delim("summaries_modules_level4.txt", delim = "\t")

summary %>%
  getSummary("modules", "level5") %>%
  write_delim("summaries_modules_level5.txt", delim = "\t")
```
