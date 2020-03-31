library("tidyverse")
library("keggR")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

load("05_MAGs/PLOTS/abund.RData")

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

## Abundance
abundance <- read_delim("05_MAGs/SUMMARY/bins_across_samples/abundance.txt", delim = "\t") %>%
  arrange(bins) %>%
  select(-bins) %>%
  rename_all(funs(str_to_lower(.)))

## Taxonomy
taxonomy <- readTax()

## Gene calls
gene_calls <- read_delim("05_MAGs/gene_calls.txt", delim = "\t") %>%
  mutate(MAG = contig %>% str_extract("[bog|heath]+_MAG_[0-9]+"))

## Gene annotation
gene_annot <- read_delim("05_MAGs/gene_calls_functions.txt", delim = "\t")

## KEGG BLAST results
kegg <- readBlast("05_MAGs/KEGG_diamond.txt")

# Reorder samples
SAMPLES <- metadata %>%
  pull(Sample) %>%
  as.vector

abundance <- abundance %>%
  select(all_of(SAMPLES))

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

# Compute mean coverage by ecosystem
# abundance.eco <- bind_cols(MAG = MAGs,
#                            Barren = abundance.sum[["Barren"]] %>%
#                              arrange(MAG) %>%
#                              pull(Mean),
#                            Heathland = abundance.sum[["Heathland"]] %>%
#                              arrange(MAG) %>%
#                              pull(Mean),
#                            Wetland = abundance.sum[["Wetland"]] %>%
#                              arrange(MAG) %>%
#                              pull(Mean))


##### HEATMAP #####

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
kegg <- kegg %>%
  assignKEGG

# Split KO table by MAG
kegg <- lapply(MAGs, function(x) {
  GENE_CALLS <- gene_calls %>%
    filter(MAG == x) %>%
    pull(gene_callers_id)

  kegg %>%
    filterKOtable(GENE_CALLS)
}) %>% set_names(MAGs)

# Summarise modules
summary <- lapply(kegg, function (MAG) {
  MAG %>%
    summariseKEGG
})

# Merge summaries
summary <- summary %>%
  mergeSummaries

# Get summaries level5
genes <- summary %>%
  getSummary("modules", "level5") %>%
  mutate_at(MAGs, function (x) ifelse(x > 0, 1, 0))

# M00567: Methanogenesis, CO2 => methane
M00567 <- selectLevel4("M00567")

pdf("05_MAGs/PLOTS/M00567-methanogenesis.pdf", width = 11.69, height = 8.27)
plotModules(M00567, SCALING = 4)
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
formatGFF(M00175_contigs) %>%
  write_delim("05_MAGs/nifH.gff", delim = "\t", col_names = F)

## M00528: Nitrification, ammonia => nitrite
M00528 <- selectLevel4("M00528")

pdf("05_MAGs/PLOTS/M00528-nitrification.pdf", width = 11.69, height = 8.27)
plotModules(M00528)
devClose()

## M00529: Denitrification, nitrate => nitrogen
M00529 <- selectLevel4("M00529")

pdf("05_MAGs/PLOTS/M00529-denitrification.pdf", width = 11.69, height = 8.27)
plotModules(M00529, SCALING = 5)
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
