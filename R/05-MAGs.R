library("tidyverse")
library("keggR")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

# Create list of MAGs
MAGs <- read_lines("05_MAGs/non_redundant_MAGs.txt")

# Read sample metadata
metadata <- readMetadata("00_METADATA/metadata.tsv") %>% 
  filter(Sample %in% SAMPLES)

# Read KEGG metadata
loadKEGG("~/KEGG")

# Read data

## Coverage
coverage <- read_delim("05_MAGs/SUMMARY/bins_across_samples/mean_coverage.txt", delim = "\t") %>% 
  arrange(bins) %>% 
  select(-bins) %>% 
  rename_all(funs(str_to_lower(.)))

# coverage.RNA <- read_delim("05_MAGs/REFINED_MAGS_RNASEQ_SUMMARY/bins_across_samples/mean_coverage.txt", delim = "\t") %>% 
#   arrange(bins) %>% 
#   select(-bins) %>% 
#   rename_all(funs(str_replace_all(., "_RNASEQ", "") %>% str_to_lower(.)))

## Taxonomy
taxonomy <- readTax()

## Gene calls (WITH ANNOTATION AFTER I RE-SUMMARISED)
gene_calls <- read_delim("05_MAGs/gene_calls.txt", delim = "\t") %>% 
  mutate(MAG = contig %>% str_extract("[bog|heath]+_MAG_[0-9]+")) %>% 
  select(gene_callers_id, MAG)

## KEGG BLAST results
blast <- readBlast("05_MAGs/gene_calls_KEGG.txt")

# Assign KOs and modules to KEGG hits
KOtable <- blast %>% 
    assignKEGG

# Split KO table
KOtable <- lapply(MAGs, function(x) {
  GENE_CALLS <- gene_calls %>%
    filter(MAG == x) %>%
    pull(gene_callers_id)

  KOtable %>%
    filterKOtable(GENE_CALLS)
}) %>% set_names(MAGs)

# Get modules
modules <- lapply(MAGs, function (MAG) {
  KOtable[[MAG]] %>%
    summariseKEGG %>% 
    getSummary("modules", "level4") %>% 
    mutate(MAG = MAG) %>% 
    mutate(count = ifelse(count > 0, 1, 0))
}) %>% bind_rows %>%
  spread(MAG, count, fill = F)

# # Genes
# genes <- lapply(MAGs, function (x) {
#   GENE_CALLS <- gene_calls %>%
#     filter(MAG == x) %>%
#     pull(gene_callers_id)
# 
#   KOtable %>%
#     filter(sequence %in% GENE_CALLS) %>%
#     summariseModules("level5", "full") %>%
#     mutate(MAG = x) %>%
#     mutate(Present = ifelse(Count > 0, 1, 0)) %>% 
#     select(-Count)
# }) %>% bind_rows %>%
#   spread(MAG, Present, fill = F)

# Reorder samples
SAMPLES <- metadata %>% 
  pull(Sample) %>% 
  as.vector

coverage <- coverage %>% 
  select(all_of(SAMPLES))

# Summarise means and standard deviations
coverage.sum <- list(All = coverage %>% 
                       summariseMeans(taxonomy),
                     Barren = coverage %>% 
                       select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
                       summariseMeans(taxonomy),
                     Heathland = coverage %>% 
                       select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
                       summariseMeans(taxonomy),
                     Wetland = coverage %>% 
                       select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>% 
                       summariseMeans(taxonomy))

coverage.sum.eco <- bind_cols(MAG = MAGs,
                              Barren = coverage.sum[["Barren"]] %>% 
                                arrange(MAG) %>% 
                                pull(Mean),
                              Heathland = coverage.sum[["Heathland"]] %>% 
                                arrange(MAG) %>% 
                                pull(Mean),
                              Wetland = coverage.sum[["Wetland"]] %>% 
                                arrange(MAG) %>% 
                                pull(Mean))

# Write mean coverage by ecosystem
write_delim(coverage.sum.eco, "05_MAGs/PLOTS/mean_coverage_ecoystem.txt", delim = "\t")


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
coverage.rich <- coverage %>% 
  computeRich

# Plot
png("05_MAGs/PLOTS/COVERAGE-richness.png", width = 1600, height = 1200, res = 150)
plotBoxplot(coverage.rich, Ecosystem)


##### ORDINATION #####

library("vegan")

# Prepare data
coverage.ord <- coverage %>% 
  sqrt

# NP-MANOVA
adonis(coverage.ord %>% t ~ Layer*Ecosystem, metadata, permutations = 9999, method = "bray")

# NMDS
coverage.mds <- coverage %>% 
  t %>% 
  metaMDS(distance = "bray")

coverage.mds.min <- coverage %>% 
  select(metadata %>% filter(Layer == "mineral") %>% pull(Sample)) %>% 
  t %>% 
  metaMDS(distance = "bray")

coverage.mds.org <- coverage %>% 
  select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  t %>% 
  metaMDS(distance = "bray")

# Plot
png("05_MAGs/PLOTS/COVERAGE-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdination(coverage.mds, metadata, "Ecosystem")

png("05_MAGs/R-PLOTS/COVERAGE-NMDS-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(coverage.mds.min, metadata %>% filter(Layer == "mineral"), "Ecosystem")

png("05_MAGs/R-PLOTS/COVERAGE-NMDS-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(coverage.mds.org, metadata %>% filter(Layer == "organic"), "Ecosystem")
metadata %>% 
  filter(Layer == "organic") %>% 
  select(all_of(FLUX.VARS)) %>% 
  envfit(coverage.mds.org, .) %>% 
  plot(col = "#b3a5cb")


##### HEATMAP KEGG MODULES #####

library("cowplot")
library("ggplotify")
library("pheatmap")

methane <- selectLevel3("Methane metabolism")
nitrogen <- selectLevel3("Nitrogen metabolism")
carbohydrate <- selectLevel3("Central carbohydrate metabolism")
carbon <- selectLevel3("Carbon fixation")

pdf("05_MAGs/PLOTS/methane.pdf", width = 6, height = 80)
plotHeatmapKEGG(methane, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/nitrogen.pdf", width = 6, height = 80)
plotHeatmapKEGG(nitrogen, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/carbohydrate.pdf", width = 6, height = 80)
plotHeatmapKEGG(carbohydrate, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/carbon.pdf", width = 6, height = 80)
plotHeatmapKEGG(carbon, c(1, 10))
devClose()


##### INVESTIGATE MODULES #####

library("cowplot")
library("ggplotify")
library("pheatmap")

M00174 <- selectLevel4("M00174")
M00175 <- selectLevel4("M00175")
M00528 <- selectLevel4("M00528")
M00529 <- selectLevel4("M00529")

pdf("05_MAGs/PLOTS/M00174.pdf", width = 6, height = 10)
plotHeatmapKEGG2(M00174, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/M00175.pdf", width = 6, height = 10)
plotHeatmapKEGG2(M00175, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/M00528.pdf", width = 6, height = 10)
plotHeatmapKEGG2(M00528, c(1, 10))
devClose()

pdf("05_MAGs/PLOTS/M00529.pdf", width = 6, height = 25)
plotHeatmapKEGG2(M00529, c(1, 10))
devClose()


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
