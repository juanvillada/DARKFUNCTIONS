# Initialize
library("tidyverse")
source("~/Data/Scripts/DARKFUNCTIONS/R/00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

# Read metadata
metadata <- readMetadata() %>% 
  filter(Sample %in% SAMPLES)

# Read data
modules <- read_delim("02_KEGG/summaries_modules_level4.txt", delim = "\t")

# Split into counts and pathway hierarchy
modules.hier <- modules %>% 
  select(-SAMPLES)

modules <- modules %>% 
  select(SAMPLES)

# Reorder samples
SAMPLES <- metadata %>% 
  pull(Sample) %>% 
  as.vector

modules <- modules %>% 
  select(SAMPLES)

# Compute number of reads per sample
nreads <- tibble(Sample = SAMPLES,
                 Reads = colSums(modules)) %>% 
  mutate(Proportion = (Reads / 2000000) * 100)

# Transform to relative abundance
modules.rel <- modules %>% 
  mutate_all(function (x) x / sum(x))

# Summarise means and standard deviations
modules.sum <- modules.rel %>% 
  summariseMeans(modules.hier)

modules.sum.barren <- modules.rel %>% 
  select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
  summariseMeans(modules.hier)

modules.sum.heath <- modules.rel %>% 
  select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
  summariseMeans(modules.hier)

modules.sum.wet <- modules.rel %>% 
  select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>% 
  summariseMeans(modules.hier)


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
modules.rich <- modules %>% 
  computeRich()

# Plot
png("02_KEGG/KEGG-richness.png", width = 1600, height = 1200, res = 150)
plotBoxplot(modules.rich, Ecosystem)
devClose()


##### HEATMAP #####

library("pheatmap")

# Prepare data
modules.map <- bind_cols(modules.hier, modules.rel)

# Keep only pathways with relative abundance > 0.1%
modules.map <- modules.map %>% 
  filter(`Level 4` %in% (modules.sum %>% filter(Mean > 0.001) %>% pull(`Level 4`) %>% as.vector))

# Reorder taxa
modules.map <- modules.map %>% 
  arrange(match(`Level 4`, modules.sum %>% pull(`Level 4`))) %>% 
  arrange(`Level 2`, `Level 3`)

# Transform to data.frame
modules.map <- data.frame(modules.map, row.names = rownames(modules.map), check.names = F)

# Plot
plotHeatmap(modules.map,
            filename = "02_KEGG/KEGG-heatmap.png",
            gaps_row = table(modules.map$`Level 2`) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
            annotation_row = modules.map %>% select(`Level 2`),
            annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                   `Cellular processes` = "#feb8a6",
                                                   `Energy metabolism` = "#82e5ed",
                                                   `Environmental information processing` = "#e4b4e2",
                                                   `Genetic information processing` = "#a0eacd",
                                                   `Metabolism` = "#ffe0e7",
                                                   `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                   `Secondary metabolism` = "#aeac9b")),
            labels_row = paste(modules.map %>% pull(`Level 3`),
                               modules.map %>% pull(`Level 4`), sep = " / ") %>% as.character)
devClose()

# Keep only "Carbohydrate and lipid metabolism","Energy metabolism" and "Nucleotide and amino acid metabolism"
modules.map.small <- modules.map %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

modules.map.small <- data.frame(modules.map.small, row.names = rownames(modules.map.small), check.names = F)

plotHeatmap(modules.map.small,
            filename = "02_KEGG/KEGG-heatmap-small.png",
            gaps_row = table(modules.map.small$`Level 2`) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
            annotation_row = modules.map.small %>% select(`Level 2`),
            annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                   `Energy metabolism` = "#82e5ed",
                                                   `Nucleotide and amino acid metabolism` = "#94b2a6")),
            labels_row = paste(modules.map.small %>% pull(`Level 3`),
                               modules.map.small %>% pull(`Level 4`), sep = " / ") %>% as.character)
devClose()


##### ORDINATION #####

library("vegan")

# Prepare data
modules.ord <- modules.rel %>% 
  sqrt

# NP-MANOVA
adonis(modules.ord %>% t ~ Layer*Ecosystem, metadata, permutations = 9999, method = "bray")

# NMDS
modules.mds <- modules.ord %>% 
  t %>% 
  metaMDS(distance = "bray")

modules.mds.min <- modules.ord %>% 
  select(metadata %>% filter(Layer == "mineral") %>% pull(Sample)) %>% 
  t %>%  
  metaMDS(distance = "bray")

modules.mds.org <- modules.ord %>% 
  select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  t %>% 
  metaMDS(distance = "bray")

# Plot
png("02_KEGG/KEGG-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdination(modules.mds, metadata, "Ecosystem")
devClose()

png("02_KEGG/KEGG-NMDS-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(modules.mds.min, metadata %>% filter(Layer == "mineral"), "Ecosystem")
devClose()

png("02_KEGG/KEGG-NMDS-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(modules.mds.org, metadata %>% filter(Layer == "organic"), "Ecosystem")
metadata %>% 
  filter(Layer == "organic") %>% 
  select(FLUX.VARS) %>% 
  envfit(modules.mds.org, .) %>% 
  plot(col = "#b3a5cb")
devClose()


##### NEGATIVE BINOMIAL ANALYSES (VEGETATION) #####

library("DESeq2")
library("pheatmap")

# Run models
modules.bin.veg <- runDeseq(modules, metadata, "Ecosystem")

# Get normalised counts
modules.bin.veg.norm <- counts(modules.bin.veg, normalized = T) %>% 
  as_tibble

# Get results
modules.bin.veg <- results(modules.bin.veg) %>% 
  as_tibble

# Plot

## Prepare data
modules.bin.veg.map <- bind_cols(modules.bin.veg.norm, modules.bin.veg, modules.hier)

## Keep only p < 0.05, reorder by log2FoldChange AND Keep only "Carbohydrate and lipid metabolism","Energy metabolism" and "Nucleotide and amino acid metabolism"
modules.bin.veg.map <- modules.bin.veg.map %>% 
  filterDeseq(0.05) %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

## Transform to data.frame
modules.bin.veg.map <- data.frame(modules.bin.veg.map, row.names = rownames(modules.bin.veg.map), check.names = F)

## Plot
plotHeatmapBin(modules.bin.veg.map, 
               filename = "02_KEGG/KEGG-heatmap-bin-ecosystem.png",
               annotation_row = modules.bin.veg.map %>% select(`Level 2`),
               annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                        Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                        `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                      `Energy metabolism` = "#82e5ed",
                                                      `Nucleotide and amino acid metabolism` = "#94b2a6")),
               labels_row = paste(modules.bin.veg.map %>% pull(`Level 3`),
                                  modules.bin.veg.map %>% pull(`Level 4`), sep = " | ") %>% as.character)
devClose()


##### NEGATIVE BINOMIAL ANALYSES (FLUX) #####

library("DESeq2")
library("pheatmap")

# Run models
modules.bin.ch <- runDeseq(modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "CH4")
modules.bin.co <- runDeseq(modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "CO2")
modules.bin.no <- runDeseq(modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "N2O")

# Get normalised counts
modules.bin.ch.norm <- counts(modules.bin.ch, normalized = T) %>% 
  as_tibble

modules.bin.co.norm <- counts(modules.bin.co, normalized = T) %>% 
  as_tibble

modules.bin.no.norm <- counts(modules.bin.no, normalized = T) %>% 
  as_tibble

# Get results
modules.bin.ch <- results(modules.bin.ch) %>% 
  as_tibble

modules.bin.co <- results(modules.bin.co) %>% 
  as_tibble

modules.bin.no <- results(modules.bin.no) %>% 
  as_tibble

# Plot

## Prepare data
modules.bin.ch.map <- bind_cols(modules.bin.ch.norm, modules.bin.ch, modules.hier)
modules.bin.co.map <- bind_cols(modules.bin.co.norm, modules.bin.co, modules.hier)
modules.bin.no.map <- bind_cols(modules.bin.no.norm, modules.bin.no, modules.hier)

## Keep only p < 0.05, reorder by log2FoldChange AND Keep only "Carbohydrate and lipid metabolism","Energy metabolism" and "Nucleotide and amino acid metabolism"
modules.bin.ch.map <- modules.bin.ch.map %>% 
  filterDeseq(0.05) %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

modules.bin.co.map <- modules.bin.co.map %>% 
  filterDeseq(0.05) %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

modules.bin.no.map <- modules.bin.no.map %>% 
  filterDeseq(0.05) %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

## Transform to data.frame
modules.bin.ch.map <- data.frame(modules.bin.ch.map, row.names = rownames(modules.bin.ch.map), check.names = F)
modules.bin.co.map <- data.frame(modules.bin.co.map, row.names = rownames(modules.bin.co.map), check.names = F)
modules.bin.no.map <- data.frame(modules.bin.no.map, row.names = rownames(modules.bin.no.map), check.names = F)

# Plot
plotHeatmapBin2(modules.bin.ch.map, "CH4",
                filename = "02_KEGG/KEGG-heatmap-bin-CH4.png",
                annotation_row = modules.bin.ch.map %>% select(`Level 2`),
                annotation_colors = list(`Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6")),
                labels_row = paste(modules.bin.ch.map %>% pull(`Level 3`),
                                   modules.bin.ch.map %>% pull(`Level 4`), sep = " | ") %>% as.character)
devClose()

plotHeatmapBin2(modules.bin.co.map, "CO2",
                filename = "02_KEGG/KEGG-heatmap-bin-CO2.png",
                annotation_row = modules.bin.co.map %>% select(`Level 2`),
                annotation_colors = list(`Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6")),
                labels_row = paste(modules.bin.co.map %>% pull(`Level 3`),
                                   modules.bin.co.map %>% pull(`Level 4`), sep = " | ") %>% as.character)
devClose()


##### OLD (WITH NUMERICAL VARIABLES) #####

##### ORDINATION

# Create list of samples with full metadata
SAMPLES.FULL <- metadata %>% 
  select(Sample, NUM.VARS) %>% 
  drop_na %>% 
  pull(Sample) %>% 
  as.vector

# NP-MANOVA
adonis(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., 
       metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), permutations = 9999, method = "bray")

# dbRDA
modules.rda <- ordiR2step(dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), 
                          dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), "forward")

png("02_KEGG/KEGG-dbRDA-layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(modules.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
devClose()

png("02_KEGG/KEGG-dbRDA-vegetation.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(modules.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
devClose()


##### NEGATIVE BINOMIAL MODELLING

library("DESeq2")

# Run models
modules.bin.ph <- runDeseq(modules %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "pH")
modules.bin.wt <- runDeseq(modules %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SoilWater")
modules.bin.om <- runDeseq(modules %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SOM")

# Get normalised counts
modules.norm <- counts(modules.bin.ph, normalized = T) %>% 
  as_tibble

# Get results
modules.bin.ph <- results(modules.bin.ph) %>% 
  as_tibble

modules.bin.wt <- results(modules.bin.wt) %>% 
  as_tibble

modules.bin.om <- results(modules.bin.om) %>% 
  as_tibble


##### HEATMAP NEGATIVE BINOMIAL MODELLING

library("pheatmap")

# Prepare data
modules.bin.ph.map <- bind_cols(modules.norm, modules.bin.ph, modules.hier)
modules.bin.wt.map <- bind_cols(modules.norm, modules.bin.wt, modules.hier)
modules.bin.om.map <- bind_cols(modules.norm, modules.bin.om, modules.hier)

# Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
modules.bin.ph.map <- filterDeseq2(modules.bin.ph.map)
modules.bin.wt.map <- filterDeseq2(modules.bin.wt.map)
modules.bin.om.map <- filterDeseq2(modules.bin.om.map)

# Transform to data.frame
modules.bin.ph.map <- data.frame(modules.bin.ph.map, row.names = rownames(modules.bin.ph.map), check.names = F)
modules.bin.wt.map <- data.frame(modules.bin.wt.map, row.names = rownames(modules.bin.wt.map), check.names = F)
modules.bin.om.map <- data.frame(modules.bin.om.map, row.names = rownames(modules.bin.om.map), check.names = F)

# Plot
plotHeatmapBin(modules.bin.ph.map, "pH",
               filename = "02_KEGG/KEGG-pH-heatmap.png",
               annotation_row = cbind(modules.bin.ph.map["Level 2"], modules.bin.ph.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                      `Cellular processes` = "#feb8a6",
                                                      `Energy metabolism` = "#82e5ed",
                                                      `Environmental information processing` = "#e4b4e2",
                                                      `Genetic information processing` = "#a0eacd",
                                                      `Metabolism` = "#ffe0e7",
                                                      `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                      `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.ph.map$`Level 3`,
                                               modules.bin.ph.map$`Level 4`, sep = " / ")))

plotHeatmapBin(modules.bin.wt.map, "SoilWater",
               filename = "02_KEGG/KEGG-SoilWater-heatmap.png",
               annotation_row = cbind(modules.bin.wt.map["Level 2"], modules.bin.wt.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                      `Cellular processes` = "#feb8a6",
                                                      `Energy metabolism` = "#82e5ed",
                                                      `Environmental information processing` = "#e4b4e2",
                                                      `Genetic information processing` = "#a0eacd",
                                                      `Metabolism` = "#ffe0e7",
                                                      `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                      `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.wt.map$`Level 3`,
                                               modules.bin.wt.map$`Level 4`, sep = " / ")))

plotHeatmapBin(modules.bin.om.map, "SOM",
               filename = "02_KEGG/KEGG-SOM-heatmap.png",
               annotation_row = cbind(modules.bin.om.map["Level 2"], modules.bin.om.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                      `Cellular processes` = "#feb8a6",
                                                      `Energy metabolism` = "#82e5ed",
                                                      `Environmental information processing` = "#e4b4e2",
                                                      `Genetic information processing` = "#a0eacd",
                                                      `Metabolism` = "#ffe0e7",
                                                      `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                      `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.om.map$`Level 3`,
                                               modules.bin.om.map$`Level 4`, sep = " / ")))
