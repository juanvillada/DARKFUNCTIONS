library("tidyverse")
library("keggR")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

load("05_MAGs/PLOTS/R.RData")

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

SAMPLES.RNA <- read_lines("SAMPLES_RNAseq.txt") %>% 
  gsub("_RNAseq", "", .)

# Create list of MAGs
MAGs <- read_lines("05_MAGs/non_redundant_MAGs.txt")

# Read sample metadata
metadata <- readMetadata("00_METADATA/metadata.tsv") %>%
  filter(Sample %in% SAMPLES)

metadata_RNA <- readMetadata("00_METADATA/metadata.tsv") %>% 
  filter(Sample %in% SAMPLES.RNA)

# Read KEGG metadata
loadKEGG("~/KEGG")

# Read data

## MAG abundance
abundance <- read_delim("05_MAGs/SUMMARY/bins_across_samples/abundance.txt", delim = "\t") %>%
  arrange(bins) %>%
  select(-bins) %>%
  rename_all(funs(str_to_lower(.)))

## Gene coverage in RNAseq
coverage <- lapply(MAGs, function(MAG) {
  read_delim(paste("05_MAGs/SUMMARY_RNASEQ/bin_by_bin/", MAG, "/" , MAG, "-gene_coverages.txt", sep = ""), delim = "\t")
}) %>% 
  bind_rows %>% 
  rename_all(funs(str_replace_all(., "_RNASEQ", "") %>% str_to_lower(.))) %>% 
  arrange(gene_callers_id)

## Taxonomy
taxonomy <- readTax()

## Gene calls
gene_calls <- read_delim("05_MAGs/gene_calls.txt", delim = "\t") %>%
  mutate(MAG = contig %>% str_extract("[bog|heath]+_MAG_[0-9]+"))

## Gene annotation
gene_annot <- read_delim("05_MAGs/gene_calls_functions.txt", delim = "\t")

## KEGG annotation
kegg <- list(BLAST = readBlast("05_MAGs/KEGG_diamond.txt", e_value = T),
             HMM = readHMM("05_MAGs/KOFAM/KOFAM_table.txt", e_value = T))

# Reorder samples
SAMPLES <- metadata %>%
  pull(Sample) %>%
  as.vector

SAMPLES.RNA <- metadata_RNA %>% 
  pull(Sample) %>% 
  as.vector

abundance <- abundance %>%
  select(all_of(SAMPLES))

coverage <- coverage %>% 
  select(all_of(SAMPLES.RNA))

# Transform to relative abundance
abundance.rel <- abundance %>%
  mutate_all(function (x) x / sum(x))

# Summarise means and standard deviations
abundance.sum <- list(All = abundance.rel %>%
                        summariseMeans(taxonomy),
                      Barren = abundance.rel %>%
                        select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>%
                        summariseMeans(taxonomy),
                      Heathland = abundance.rel %>%
                        select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>%
                        summariseMeans(taxonomy),
                      Wetland = abundance.rel %>%
                        select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>%
                        summariseMeans(taxonomy))

coverage.sum <- list(All = coverage %>% 
                       summariseMeans(gene_calls),
                     Barren = coverage %>% 
                       select(metadata_RNA %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
                       summariseMeans(gene_calls),
                     Heathland = coverage %>% 
                       select(metadata_RNA %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
                       summariseMeans(gene_calls))


##### HEATMAP MAG ABUNDANCE #####

library("pheatmap")

# Prepare data
abundance.map <- bind_cols(taxonomy, abundance.rel)

# Reorder taxa
abundance.map <- abundance.map %>%
  arrange(Domain, Phylum, MAG)

# Plot
plotHeatmap(abundance.map,
            filename = "05_MAGs/PLOTS/ABUNDANCE-heatmap.png",
            gaps_row = table(abundance.map$Domain) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
            annotation_row = data.frame(abundance.map %>% select(Phylum), row.names = rownames(abundance.map)),
            annotation_colors = list(),
            labels_row = abundance.map %>% pull(MAG) %>% as.character)
devClose()


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
abundance.rich <- abundance %>%
  computeRich

# Plot
png("05_MAGs/PLOTS/ABUNDANCE-richness.png", width = 1600, height = 1200, res = 150)
plotBoxplot(abundance.rich, Ecosystem)
devClose()


##### ORDINATION #####

library("vegan")

# Prepare data
abundance.ord <- abundance %>%
  sqrt

# NP-MANOVA
adonis(abundance.ord %>% t ~ Layer*Ecosystem, metadata, permutations = 9999, method = "bray")

# NMDS
abundance.mds <- abundance %>%
  t %>%
  metaMDS(distance = "bray")

abundance.mds.min <- abundance %>%
  select(metadata %>% filter(Layer == "mineral") %>% pull(Sample)) %>%
  t %>%
  metaMDS(distance = "bray")

abundance.mds.org <- abundance %>%
  select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>%
  t %>%
  metaMDS(distance = "bray")

# Plot
png("05_MAGs/PLOTS/ABUNDANCE-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdination(abundance.mds, metadata, "Ecosystem")
devClose()

png("05_MAGs/PLOTS/ABUNDANCE-NMDS-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(abundance.mds.min, metadata %>% filter(Layer == "mineral"), "Ecosystem")
devClose()

png("05_MAGs/PLOTS/ABUNDANCE-NMDS-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(abundance.mds.org, metadata %>% filter(Layer == "organic"), "Ecosystem")
metadata %>%
  filter(Layer == "organic") %>%
  select(all_of(FLUX.VARS)) %>%
  envfit(abundance.mds.org, .) %>%
  plot(col = "#b3a5cb")
devClose()


##### KEGG MODULES #####

library("cowplot")
library("ggplotify")
library("pheatmap")

# Assign KOs and modules to KEGG hits
kegg[["BLAST"]] <- kegg[["BLAST"]] %>% 
  assignKEGG

kegg[["HMM"]] <- kegg[["HMM"]] %>% 
  assignKEGG

# Split KO table by MAG
kegg <- lapply(kegg, function (x) {
  lapply(MAGs, function(y) {
    GENE_CALLS <- gene_calls %>%
      filter(MAG == y) %>%
      pull(gene_callers_id)
    
    x %>%
      filterKOtable(GENE_CALLS)
  }) %>% set_names(MAGs)
})

# Summarise modules
summary <- lapply(kegg, function(x) {
  lapply(x, function (MAG) {
    MAG %>%
      summariseKEGG
  })
})

# Merge summaries
summary <- lapply(summary, function(x) {
  x %>% 
    mergeSummaries
})

# Get summaries level5
genes <- lapply(summary, function(x) {
  x %>% 
    getSummary("modules", "level5") 
})

genes <- genes %>% 
  bind_rows %>% 
  group_by(level1, level2, level3, level4, KO, gene) %>% 
  summarise_all(funs(sum(., na.rm = T))) %>% 
  ungroup %>%
  mutate_at(MAGs, function (x) ifelse(x > 0, 1, 0))

# M00567: Methanogenesis, CO2 => methane
M00567 <- selectLevel4("M00567")

pdf("05_MAGs/PLOTS/M00567-methanogenesis.pdf", width = 11.69, height = 8.27)
plotModules(M00567, SCALING = 2.8)
devClose()

# M00174: Methane oxidation, methanotroph, methane => formaldehyde
M00174 <- selectLevel4("M00174")

pdf("05_MAGs/PLOTS/M00174-methane-oxidation.pdf", width = 11.69, height = 8.27)
plotModules(M00174)
devClose()

## Investigate contigs
extractGene("bog_MAG_0057", "K16157")
extractContig("bog_MAG_0057_000000000172")

# M00175: Nitrogen fixation, nitrogen => ammonia
M00175 <- selectLevel4("M00175")

## Add some other nif genes
M00175 <- M00175 %>%
  addToLevel4(c("K02584", "K02585", "K02587", "K02589", "K02590", "K02592", "K02593", "K02594", "K02595", "K02596", "K02597", "K04488", "K15790"))

pdf("05_MAGs/PLOTS/M00175-nitrogen-fixation.pdf", width = 11.69, height = 8.27)
plotModules(M00175)
devClose()

## Investigate contigs
M00175_genes <- lapply(M00175$MAGs, function(x) {
  extractGene(x, c("K02584", "K02585", "K02586", "K02587", "K02588", "K02589", "K02590", "K02591", "K02592", "K02593", "K02594", "K02595", "K02596", "K02597", "K04488", "K15790"))
}) %>%set_names(M00175$MAGs)

viewList(M00175_genes)

M00175_contigs <- lapply(M00175$MAGs, function(x) {
  M00175_genes[[x]][["gene_calls"]]
}) %>%
  bind_rows %>%
  pull(contig) %>%
  unique

M00175_contigs <- lapply(M00175_contigs, function(x) {
  extractContig(x)
}) %>%set_names(M00175_contigs)

viewList(M00175_contigs)

## Create GFF file
formatGFF(M00175_contigs, "KEGG") %>%
  write_delim("05_MAGs/TEST_NIF/nifH-KEGG.gff", delim = "\t", col_names = F)

formatGFF(M00175_contigs, "KOFAM") %>%
  write_delim("05_MAGs/TEST_NIF/nifH-KOFAM.gff", delim = "\t", col_names = F)

## M00528: Nitrification, ammonia => nitrite
M00528 <- selectLevel4("M00528")

pdf("05_MAGs/PLOTS/M00528-nitrification.pdf", width = 11.69, height = 8.27)
plotModules(M00528)
devClose()

## M00529: Denitrification, nitrate => nitrogen
M00529 <- selectLevel4("M00529")

pdf("05_MAGs/PLOTS/M00529-denitrification.pdf", width = 11.69, height = 8.27)
plotModules(M00529, SCALING = 4)
devClose()


##### COVERAGE RNA ######

library("pheatmap")

# Select top 100 most transcribed genes
coverage.top100 <- coverage.sum[["All"]] %>% 
  filter(source != "Ribosomal_RNAs") %>% 
  slice(1:100) %>% 
  pull(gene_callers_id)

coverage.top100 <- bind_cols(gene_calls, coverage) %>%
  filter(gene_callers_id %in% coverage.top100) %>% 
  select(gene_callers_id, MAG, all_of(SAMPLES.RNA)) %>% 
  left_join(gene_annot %>% filter(source == "COG_CATEGORY")) %>% 
  arrange(MAG)

# Heatmap
plotHeatmapCoverage <- function(input, ...) {
  ORDER <- metadata_RNA %>% 
    arrange(Ecosystem, Layer) %>% 
    pull(Sample)
  
  ANNOTATION_COL <- data.frame(Ecosystem = metadata_RNA %>% select(Ecosystem), row.names = metadata_RNA %>% pull(Sample))
  
  pheatmap(input %>% select(all_of(ORDER)) %>% sqrt, 
           color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100), 
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, 
           annotation_col = ANNOTATION_COL,
           # annotation_row = DATA %>% select(MAG),
           annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                    Layer = c(mineral = "#b7d8ff", organic = "#98c699")),
           labels_row = paste(input %>% pull(gene_callers_id), input %>% pull(MAG), input %>% pull(`function`), sep = " | ") %>% as.character, ...)
}

plotHeatmapCoverage(coverage.top100,
                    filename = "05_MAGs/PLOTS/GENE-COVERAGE-TOP100.png")
devClose()

# nif <- gene_calls_annot %>% 
#   filter(source == "KEGG") %>% 
#   filter(str_detect(`function`, "nif")) %>% 
#   pull(gene_callers_id)
# 
# nif <- bind_cols(gene_calls, coverage) %>% 
#   filter(gene_callers_id %in% nif)


##### CUSTOM MODULES #####

fructose <- c("K00847", "K00844")
galactose <- c("K01785", "K00849", "K00965", "K01784")

fructose <- genes %>%
  filter(KO %in% fructose) %>%
  select(-starts_with("Level")) %>%
  unique %>%
  select(-KO) %>%
  column_to_rownames("gene")

galactose <- genes %>%
  filter(KO %in% galactose) %>%
  select(-starts_with("Level")) %>%
  unique %>%
  select(-KO) %>%
  column_to_rownames("gene")

pheatmap(fructose, cluster_rows = F, cluster_cols = F)
pheatmap(galactose)
