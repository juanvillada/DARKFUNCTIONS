library("dplyr")
source("00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Set working directory
setwd("~/Helsinki/analyses/")

# Read metadata
metadata <- read.table("00_METADATA/metadata_clean.tsv", header = T, sep = "\t")

# Create list of samples
SAMPLES <- scan("SAMPLES.txt", character(), quote = "") 

# Subset metadata
metadata <- metadataSub(SAMPLES)

# Read data
modules <- read.delim("02_KEGG/summaries_modules_level4.txt", header = T, sep = "\t")

# Split into counts and pathway hierarchy
modules.hier <- modules %>% 
  select(-SAMPLES)

modules <- modules %>% 
  select(SAMPLES)

# Reorder samples
SAMPLES <- metadata %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

modules <- modules[SAMPLES]

# Transform to relative abundance
modules.rel <- sweep(modules, 2, colSums(modules), "/")


##### HEATMAP #####

library("pheatmap")

# Prepare data
modules.map <- cbind(modules.hier, modules.rel)

# Keep only pathways with relative abundance > 0.1%
modules.map <- modules.map %>% 
  filter(apply(modules.map[SAMPLES], 1, mean) > 0.001)

# Reorder taxa
modules.map <- modules.map %>% 
  arrange(desc(apply(modules.map[SAMPLES], 1, mean))) %>% 
  arrange(Level.2, Level.3)

# Fix a weird bug
modules.map <- data.frame(modules.map, row.names = rownames(modules.map), check.names = F)

# Plot
plotHeatmap(modules.map,
            filename = "02_KEGG/KEGG-heatmap.png",
            annotation_row = modules.map["Level.2"],
            annotation_colors = list(Habitat = c(`Betula nana` = "#dfc3f8", Bog = "#beefc1",
                                                 Empetrum = "#eca6c1", Heath = "#61b7d9",
                                                 Meadow = "#f9b99f", Shrub = "#d8deff",
                                                 Undefined = "#b3a5cb"),
                                     Layer = c(Mineral = "#b7d8ff", Organic = "#98c699"),
                                     Level.2 = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                 `Cellular processes` = "#feb8a6",
                                                 `Energy metabolism` = "#82e5ed",
                                                 `Environmental information processing` = "#e4b4e2",
                                                 `Genetic information processing` = "#a0eacd",
                                                 `Metabolism` = "#ffe0e7",
                                                 `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                 `Secondary metabolism` = "#aeac9b")),
            labels_row = as.character(paste(modules.map$Level.3,
                                            modules.map$Level.4, sep = " / ")))

dev.off()


##### ORDINATION #####

library("vegan")

# Create list of samples with full metadata
SAMPLES.FULL <- metadata %>% 
  na.omit %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

# Create list of numeric variables
NUM.VARS = metadata %>% 
  select(c(5:18)) %>% 
  names

# Transform numeric variables to log
metadata[NUM.VARS] <- log(metadata[NUM.VARS])

# NP-MANOVA
modules.pmn <- adonis(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), permutations = 9999, method = "bray")

# NMDS
modules.ord <- metaMDS(modules.rel %>% t %>% sqrt, distance = "bray")

# NMDS
modules.ord.full <- metaMDS(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt, distance = "bray")

# dbRDA
modules.rda <- ordiR2step(dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), 
                          dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), "forward")

# Plot
png("02_KEGG/KEGG-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(modules.ord, metadata)
dev.off()

png("02_KEGG/KEGG-NMDS_full.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(modules.ord.full, metadataSub(SAMPLES.FULL))
dev.off()

png("02_KEGG/KEGG-dbRDA_habitat.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(modules.rda, metadataSub(SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
dev.off()

png("02_KEGG/KEGG-dbRDA_layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(modules.rda, metadataSub(SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
dev.off()


##### NEGATIVE BINOMIAL MODELLING #####

library("DESeq2")

# Run models
modules.bin.ph <- DESeqDataSetFromMatrix(modules %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ pH) %>%
  DESeq(fitType = "local", quiet = T)

modules.bin.den <- DESeqDataSetFromMatrix(modules %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ Density) %>%
  DESeq(fitType = "local", quiet = T)

modules.bin.dwr <- DESeqDataSetFromMatrix(modules %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ DryWetRatio) %>%
  DESeq(fitType = "local", quiet = T)

modules.bin.som <- DESeqDataSetFromMatrix(modules %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ SOM) %>%
  DESeq(fitType = "local", quiet = T)

# Get normalised counts
modules.norm <- counts(modules.bin.ph, normalized = T) %>% 
  as.data.frame

# Get results
modules.bin.ph <- results(modules.bin.ph) %>% 
  as.data.frame

modules.bin.den <- results(modules.bin.den) %>% 
  as.data.frame

modules.bin.dwr <- results(modules.bin.dwr) %>% 
  as.data.frame

modules.bin.som <- results(modules.bin.som) %>% 
  as.data.frame


##### HEATMAP NEGATIVE BINOMIAL MODELLING #####

library("pheatmap")

# Prepare data
modules.bin.ph.map  <- cbind(modules.norm, modules.bin.ph,  modules.hier)
modules.bin.den.map <- cbind(modules.norm, modules.bin.den, modules.hier)
modules.bin.dwr.map <- cbind(modules.norm, modules.bin.dwr, modules.hier)
modules.bin.som.map <- cbind(modules.norm, modules.bin.som, modules.hier)

# Keep only p < 0.05 and reorder by log2FoldChange
modules.bin.ph.map <- modules.bin.ph.map %>% 
  subset(pvalue < 0.05) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change, Level.2)

modules.bin.den.map <- modules.bin.den.map %>% 
  subset(pvalue < 0.05) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change, Level.2)

modules.bin.dwr.map <- modules.bin.dwr.map %>% 
  subset(pvalue < 0.05) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change, Level.2)

modules.bin.som.map <- modules.bin.som.map %>% 
  subset(pvalue < 0.05 ) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change, Level.2)

# Fix a weird bug
modules.bin.ph.map  <- data.frame(modules.bin.ph.map,  row.names = rownames(modules.bin.ph.map),  check.names = F)
modules.bin.den.map <- data.frame(modules.bin.den.map, row.names = rownames(modules.bin.den.map), check.names = F)
modules.bin.dwr.map <- data.frame(modules.bin.dwr.map, row.names = rownames(modules.bin.dwr.map), check.names = F)
modules.bin.som.map <- data.frame(modules.bin.som.map, row.names = rownames(modules.bin.som.map), check.names = F)

# Plot
plotHeatmapBin(modules.bin.ph.map, "pH",
               filename = "02_KEGG/KEGG-pH-heatmap.png",
               annotation_row = cbind(modules.bin.ph.map["Level.2"], modules.bin.ph.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        Level.2 = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                    `Cellular processes` = "#feb8a6",
                                                    `Energy metabolism` = "#82e5ed",
                                                    `Environmental information processing` = "#e4b4e2",
                                                    `Genetic information processing` = "#a0eacd",
                                                    `Metabolism` = "#ffe0e7",
                                                    `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                    `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.ph.map$Level.3,
                                               modules.bin.ph.map$Level.4, sep = " / ")))

plotHeatmapBin(modules.bin.den.map, "Density",
               filename = "02_KEGG/KEGG-den-heatmap.png",
               annotation_row = cbind(modules.bin.den.map["Level.2"], modules.bin.den.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        Level.2 = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                    `Cellular processes` = "#feb8a6",
                                                    `Energy metabolism` = "#82e5ed",
                                                    `Environmental information processing` = "#e4b4e2",
                                                    `Genetic information processing` = "#a0eacd",
                                                    `Metabolism` = "#ffe0e7",
                                                    `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                    `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.den.map$Level.3,
                                               modules.bin.den.map$Level.4, sep = " / ")))

plotHeatmapBin(modules.bin.dwr.map, "DryWetRatio",
               filename = "02_KEGG/KEGG-dwr-heatmap.png",
               annotation_row = cbind(modules.bin.dwr.map["Level.2"], modules.bin.dwr.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        Level.2 = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                    `Cellular processes` = "#feb8a6",
                                                    `Energy metabolism` = "#82e5ed",
                                                    `Environmental information processing` = "#e4b4e2",
                                                    `Genetic information processing` = "#a0eacd",
                                                    `Metabolism` = "#ffe0e7",
                                                    `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                    `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.dwr.map$Level.3,
                                               modules.bin.dwr.map$Level.4, sep = " / ")))

plotHeatmapBin(modules.bin.som.map, "SOM",
               filename = "02_KEGG/KEGG-SOM-heatmap.png",
               annotation_row = cbind(modules.bin.som.map["Level.2"], modules.bin.som.map["Change"]),
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                        Level.2 = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                    `Cellular processes` = "#feb8a6",
                                                    `Energy metabolism` = "#82e5ed",
                                                    `Environmental information processing` = "#e4b4e2",
                                                    `Genetic information processing` = "#a0eacd",
                                                    `Metabolism` = "#ffe0e7",
                                                    `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                    `Secondary metabolism` = "#aeac9b")),
               labels_row = as.character(paste(modules.bin.som.map$Level.3,
                                               modules.bin.som.map$Level.4, sep = " / ")))

dev.off()

