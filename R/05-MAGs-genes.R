library("tidyverse")
library("keggR")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

load("05_MAGs/PLOTS/genes.RData")

# Create list of samples
SAMPLES <- read_lines("SAMPLES_RNAseq.txt") %>% 
  gsub("_RNAseq", "", .)

# Create list of MAGs
MAGs <- read_lines("05_MAGs/non_redundant_MAGs.txt")

# Read sample metadata
metadata <- readMetadata("00_METADATA/metadata.tsv") %>% 
  filter(Sample %in% SAMPLES)

# Read data

## Coverage
coverage <- lapply(MAGs, function(MAG) {
  read_delim(paste("05_MAGs/SUMMARY_RNASEQ/bin_by_bin/", MAG, "/" , MAG, "-gene_coverages.txt", sep = ""), delim = "\t")
}) %>% 
  bind_rows %>% 
  rename_all(funs(str_replace_all(., "_RNASEQ", "") %>% str_to_lower(.))) %>% 
  arrange(gene_callers_id)

## Gene calls
gene_calls <- read_delim("05_MAGs/gene_calls.txt", delim = "\t") %>% 
  mutate(MAG = contig %>% str_extract("[bog|heath]+_MAG_[0-9]+")) %>% 
  select(gene_callers_id, contig, MAG, source) %>% 
  arrange(gene_callers_id)

## Gene annotation
gene_calls_annot <- read_delim("05_MAGs/gene_calls_functions.txt", delim = "\t")

# Coverage

# Reorder samples
SAMPLES <- metadata %>% 
  pull(Sample) %>% 
  as.vector

coverage <- coverage %>% 
  select(all_of(SAMPLES))

# Summarise means and standard deviations
coverage.sum <- list(All = coverage %>% 
                       summariseMeans(gene_calls),
                     Barren = coverage %>% 
                       select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
                       summariseMeans(gene_calls),
                     Heathland = coverage %>% 
                       select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
                       summariseMeans(gene_calls))


top100_genes <- coverage.sum[["All"]] %>% 
  filter(source != "Ribosomal_RNAs") %>% 
  slice(1:100) %>% 
  pull(gene_callers_id)

top100_genes <- bind_cols(gene_calls, coverage) %>%
  filter(gene_callers_id %in% top100_genes) %>% 
  select(-source) %>% 
  left_join(gene_calls_annot %>% filter(source != "COG_CATEGORY")) %>% 
  arrange(MAG)

top100_genes <- data.frame(top100_genes, row.names = rownames(top100_genes))

library("pheatmap")

f2 <- function(DATA, ...) {
  ORDER <- metadata %>% arrange(Ecosystem, Layer) %>% pull(Sample) %>% 
    as.vector
  
  ANNOTATION_COL <- data.frame(Layer = metadata %>% select(Layer), 
                               Ecosystem = metadata %>% select(Ecosystem), 
                               row.names = metadata %>% pull(Sample))
  
  pheatmap(DATA %>% select(all_of(ORDER)) %>% sqrt, 
           color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100), 
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, 
           annotation_col = ANNOTATION_COL,
           # annotation_row = DATA %>% select(MAG),
           annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                    Layer = c(mineral = "#b7d8ff", organic = "#98c699")),
           labels_row = paste(DATA %>% pull(gene_callers_id), DATA %>% pull(MAG), DATA %>% pull(function.), sep = " | ") %>% as.character, ...)
}

f2(top100_genes, filename = "05_MAGs/PLOTS/GENE-COVERAGE-TOP100.png")
devClose()

nif <- gene_calls_annot %>% 
  filter(source == "KEGG") %>% 
  filter(str_detect(`function`, "nif")) %>% 
  pull(gene_callers_id)

nif <- bind_cols(gene_calls, coverage) %>% 
  filter(gene_callers_id %in% nif)
