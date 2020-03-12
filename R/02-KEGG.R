library("tidyverse")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

# Read metadata
metadata <- readMetadata("00_METADATA/metadata.tsv") %>% 
  filter(Sample %in% SAMPLES)

# Read data
modules <- read_delim("02_KEGG/summaries_modules_level4.txt", delim = "\t")

# Split into counts and pathway hierarchy
modules.hier <- modules %>% 
  select(-all_of(SAMPLES))

modules <- modules %>% 
  select(all_of(SAMPLES))

# Reorder samples
SAMPLES <- metadata %>% 
  pull(Sample) %>% 
  as.vector

modules <- modules %>% 
  select(all_of(SAMPLES))

# Compute number of reads per sample
nreads <- tibble(Sample = SAMPLES,
                 Reads = colSums(modules)) %>% 
  mutate(Proportion = (Reads / 2000000) * 100)

# Transform to relative abundance
modules.rel <- modules %>% 
  mutate_all(function (x) x / sum(x))

# Summarise means and standard deviations
modules.sum <- list(All = modules.rel %>% 
                      summariseMeans(modules.hier),
                    Barren = modules.rel %>% 
                      select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
                      summariseMeans(modules.hier),
                    Heathland = modules.rel %>% 
                      select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
                      summariseMeans(modules.hier),
                    Wetland = modules.rel %>% 
                      select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>% 
                      summariseMeans(modules.hier))


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
modules.rich <- modules %>% 
  computeRich

# Plot
png("02_KEGG/KEGG-richness.png", width = 1600, height = 1200, res = 150)
plotBoxplot(modules.rich, Ecosystem)


##### HEATMAP #####

library("pheatmap")

# Prepare data
modules.map <- bind_cols(modules.hier, modules.rel)

# Keep only pathways with relative abundance > 0.1%
modules.map <- modules.map %>% 
  filter(`Level 4` %in% (modules.sum[["All"]] %>% filter(Mean > 0.001) %>% pull(`Level 4`) %>% as.vector))

# Reorder taxa
modules.map <- modules.map %>% 
  arrange(match(`Level 4`, modules.sum[["All"]] %>% pull(`Level 4`))) %>% 
  arrange(`Level 2`, `Level 3`)

# Transform to data.frame
modules.map <- data.frame(modules.map, row.names = rownames(modules.map), check.names = F)

# Plot
plotHeatmap(modules.map,
            # filename = "02_KEGG/KEGG-heatmap.png",
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
            # filename = "02_KEGG/KEGG-heatmap-small.png",
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

png("02_KEGG/KEGG-NMDS-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(modules.mds.min, metadata %>% filter(Layer == "mineral"), "Ecosystem")

png("02_KEGG/KEGG-NMDS-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(modules.mds.org, metadata %>% filter(Layer == "organic"), "Ecosystem")
metadata %>% 
  filter(Layer == "organic") %>% 
  select(all_of(FLUX.VARS)) %>% 
  envfit(modules.mds.org, .) %>% 
  plot(col = "#b3a5cb")


##### NEGATIVE BINOMIAL ANALYSES (VEGETATION) #####

library("DESeq2")
library("pheatmap")

# Run models
modules.bin.veg <- modules %>% 
  runDeseq(metadata, "Ecosystem")

# Get normalised counts
modules.bin.veg.norm <- modules.bin.veg %>% 
  counts(normalized = T) %>% 
  as_tibble

# Get results
modules.bin.veg <- modules.bin.veg %>% 
  results %>% 
  as_tibble

# Plot

## Prepare data
modules.bin.veg.map <- bind_cols(modules.bin.veg.norm, modules.bin.veg, modules.hier)

## Keep only p < 0.05, reorder by log2FoldChange and keep only "Carbohydrate and lipid metabolism","Energy metabolism" and "Nucleotide and amino acid metabolism"
modules.bin.veg.map <- modules.bin.veg.map %>% 
  filterDeseq(0.05) %>% 
  filter(`Level 2` %in% c("Carbohydrate and lipid metabolism", "Energy metabolism", "Nucleotide and amino acid metabolism"))

## Transform to data.frame
modules.bin.veg.map <- data.frame(modules.bin.veg.map, row.names = rownames(modules.bin.veg.map), check.names = F)

## Plot
plotHeatmapBin(modules.bin.veg.map, 
               # filename = "02_KEGG/KEGG-heatmap-bin-ecosystem.png",
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
modules.bin.ch <- modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  runDeseq(metadata %>% filter(Layer == "organic"), "CH4")

modules.bin.co <- modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  runDeseq(metadata %>% filter(Layer == "organic"), "CO2")

modules.bin.no <- modules %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  runDeseq(metadata %>% filter(Layer == "organic"), "N2O")

# Get normalised counts
modules.bin.ch.norm <- modules.bin.ch %>% 
  counts(normalized = T) %>% 
  as_tibble

modules.bin.co.norm <- modules.bin.co %>% 
  counts(normalized = T) %>% 
  as_tibble

modules.bin.no.norm <- modules.bin.no %>% 
  counts(normalized = T) %>% 
  as_tibble

# Get results
modules.bin.ch <- modules.bin.ch %>% 
  results %>% 
  as_tibble

modules.bin.co <- modules.bin.co %>% 
  results %>% 
  as_tibble

modules.bin.no <- modules.bin.no %>% 
  results %>% 
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
                # filename = "02_KEGG/KEGG-heatmap-bin-CH4.png",
                annotation_row = modules.bin.ch.map %>% select(`Level 2`),
                annotation_colors = list(`Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6")),
                labels_row = paste(modules.bin.ch.map %>% pull(`Level 3`),
                                   modules.bin.ch.map %>% pull(`Level 4`), sep = " | ") %>% as.character)
devClose()

plotHeatmapBin2(modules.bin.co.map, "CO2",
                # filename = "02_KEGG/KEGG-heatmap-bin-CO2.png",
                annotation_row = modules.bin.co.map %>% select(`Level 2`),
                annotation_colors = list(`Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6")),
                labels_row = paste(modules.bin.co.map %>% pull(`Level 3`),
                                   modules.bin.co.map %>% pull(`Level 4`), sep = " | ") %>% as.character)
devClose()
