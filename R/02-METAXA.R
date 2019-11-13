library("dplyr")
library("reshape2")
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
phylum <- read.table("02_METAXA/summary_level_2.txt", header = T, sep = "\t")
genus  <- read.table("02_METAXA/summary_level_6.txt", header = T, sep = "\t")

# Split taxonomy
phylum <- data.frame(phylum[SAMPLES],  colsplit(phylum$Taxa, pattern = "\\;", names = c("Domain", "Phylum")))
genus  <- data.frame(genus[SAMPLES],   colsplit(genus$Taxa,  pattern = "\\;", names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus")))

# Keep only bacteria and archaea
phylum <- phylum %>% 
  subset(Domain %in% c("Bacteria", "Archaea"))

genus <- genus %>% 
  subset(Domain %in% c("Bacteria", "Archaea"))

# Split into counts and taxonomic hierarchy
phylum.hier <- phylum %>% 
  select(-SAMPLES)

genus.hier <- genus %>% 
  select(-SAMPLES)

phylum <- phylum %>% 
  select(SAMPLES)

genus <- genus %>% 
  select(SAMPLES)

# Reorder samples
SAMPLES <- metadata %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

phylum <- phylum[SAMPLES]
genus  <- genus [SAMPLES]

# Transform to relative abundance
phylum.rel <- sweep(phylum, 2, colSums(phylum), "/")
genus.rel  <- sweep(genus,  2, colSums(genus),  "/")

# Summarise means and standard deviations
phylum.sum <- summariseMeans(phylum.rel, phylum.hier)
genus.sum  <- summariseMeans(genus.rel,  genus.hier)


##### BARPLOT #####

library("ggplot2")
library("forcats")

# Prepare data
phylum.bar <- cbind(Phylum = phylum.hier$Phylum, phylum.rel)
genus.bar  <- cbind(Genus  = genus.hier $Genus,  genus.rel)

# Reorder taxa by abundance
phylum.bar <- phylum.bar %>% 
  arrange(desc(apply(phylum.bar[SAMPLES], 1, mean)))

genus.bar <- genus.bar %>% 
  arrange(desc(apply(genus.bar[SAMPLES], 1, mean)))

# For the genus level, remove unclassified and keep only top 30 genera
genus.bar <- genus.bar %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus)) %>%
  head(30)

# Merge with metadata
phylum.bar <- metadataMerge(phylum.bar)
genus.bar  <- metadataMerge(genus.bar)

# Plot
png("02_METAXA/METAXA-PHYLUM-barplot.png", width = 1600, height = 1200, res = 150)
plot(plotBarplot(phylum.bar)) +
  scale_fill_manual(values = c("#ffbb9b", "#47e0ff", "#fbc387", "#64cbff", "#f0c683", "#6fb1ec", "#cfdd90", "#ce9cde", "#b9f5b5", "#faa0ca",
                               "#7ff9e1", "#ffaba1", "#59f0fc", "#d1a377", "#3ab9dc", "#ffddae", "#46d0c3", "#ffbdce", "#74dcb7", "#e9ccff",
                               "#91ba76", "#c3a1bf", "#e2ffb9", "#97acd4", "#73bc87", "#f0e3ff", "#68bb9e", "#d89d9a", "#c1ffdd", "#bcaa77",
                               "#bae8ff", "#fff4c2", "#75b8ae", "#eeffe2", "#91b3a0", "#c5fff0", "#a7af8c", "#b4e1e2"))
dev.off()

png("02_METAXA/METAXA-GENUS-barplot.png", width = 1600, height = 1200, res = 150)
plot(plotBarplot(genus.bar)) +
  scale_fill_manual(values = c("#e2f0ff", "#daa883", "#89caf7", "#e6cf99", "#7fb0e1", "#f0f6bc", "#c1baf1", "#acbf89", "#d7a7d4", "#c6ffd7",
                               "#d39cb4", "#90eade", "#db9c8d", "#71d4de", "#ffc4a6", "#55bac7", "#ffcfb7", "#a7f0ff", "#cea387", "#d0ffff",
                               "#bbab7c", "#d8dbff", "#88c29d", "#ffd9f7", "#8dd8bd", "#c7dfff", "#95b396", "#eefff4", "#bba7a3", "#95b1b2"))
dev.off()


##### HEATMAP #####

library("pheatmap")

# Prepare data
phylum.map <- cbind(phylum.hier, phylum.rel)
genus.map  <- cbind(genus.hier,  genus.rel)

# For the phylum level, remove unclassified
phylum.map <- phylum.map %>% 
  subset(!grepl("Unclassified", Phylum))

# For the genus level, remove unclassified and keep only top 50 genera
genus.map <- genus.map %>% 
  arrange(desc(apply(genus.map[SAMPLES], 1, mean))) %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus)) %>%
  head(50) %>% 
  droplevels.data.frame

# Reorder taxa
phylum.map <- phylum.map %>% 
  arrange(Domain)

genus.map  <- genus.map %>% 
  arrange(Phylum, Genus)

# Fix a weird bug
phylum.map <- data.frame(phylum.map, row.names = rownames(phylum.map), check.names = F)
genus.map  <- data.frame(genus.map,  row.names = rownames(genus.map), check.names = F)

# Plot
plotHeatmap(phylum.map,
            filename = "02_METAXA/METAXA-PHYLUM-heatmap.png",
            gaps_row = table(phylum.map$Domain) %>% as.vector %>% cumsum,
            annotation_row = phylum.map["Domain"],
            annotation_colors = list(Habitat = c(`Betula nana` = "#dfc3f8", Bog = "#beefc1",
                                                 Empetrum = "#eca6c1", Heath = "#61b7d9",
                                                 Meadow = "#f9b99f", Shrub = "#d8deff",
                                                 Undefined = "#b3a5cb"),
                                     Layer = c(Mineral = "#b7d8ff", Organic = "#98c699"),
                                     Domain = c(Archaea = "#8befff", Bacteria = "#acaf79")),
            labels_row = as.character(phylum.map$Phylum))

plotHeatmap(genus.map,
            filename = "02_METAXA/METAXA-GENUS-heatmap.png",
            gaps_row = table(genus.map$Phylum) %>% as.vector %>% cumsum,
            annotation_row = genus.map["Phylum"],
            annotation_colors = list(Habitat = c(`Betula nana` = "#dfc3f8", Bog = "#beefc1",
                                                 Empetrum = "#eca6c1", Heath = "#61b7d9",
                                                 Meadow = "#f9b99f", Shrub = "#d8deff",
                                                 Undefined = "#b3a5cb"),
                                     Layer = c(Mineral = "#b7d8ff", Organic = "#98c699"),
                                     Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                Verrucomicrobia = "#93b48d")),
            labels_row = as.character(genus.map$Genus))

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
genus.pmn <- adonis(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), permutations = 9999, method = "bray")

# NMDS
genus.ord <- metaMDS(genus.rel %>% t %>% sqrt, distance = "bray")

# NMDS
genus.ord.full <- metaMDS(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt, distance = "bray")

# dbRDA
genus.rda <- ordiR2step(dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), 
                        dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadataSub(SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), "forward")

# Plot
png("02_METAXA/METAXA-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(genus.ord, metadata)
dev.off()

png("02_METAXA/METAXA-NMDS_full.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(genus.ord.full, metadataSub(SAMPLES.FULL))
dev.off()

png("02_METAXA/METAXA-dbRDA_habitat.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(genus.rda, metadataSub(SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
dev.off()

png("02_METAXA/METAXA-dbRDA_layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(genus.rda, metadataSub(SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
dev.off()


##### NEGATIVE BINOMIAL MODELLING #####

library("DESeq2")

# Run models
genus.bin.ph <- DESeqDataSetFromMatrix(genus %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ pH) %>%
  DESeq(fitType = "local", quiet = T)

genus.bin.den <- DESeqDataSetFromMatrix(genus %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ Density) %>%
  DESeq(fitType = "local", quiet = T)

genus.bin.dwr <- DESeqDataSetFromMatrix(genus %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ DryWetRatio) %>%
  DESeq(fitType = "local", quiet = T)

genus.bin.som <- DESeqDataSetFromMatrix(genus %>% select(SAMPLES.FULL), metadataSub(SAMPLES.FULL), ~ SOM) %>%
  DESeq(fitType = "local", quiet = T)

# Get normalised counts
genus.norm <- counts(genus.bin.ph, normalized = T) %>% 
    as.data.frame

# Get results
genus.bin.ph <- results(genus.bin.ph) %>% 
  as.data.frame

genus.bin.den <- results(genus.bin.den) %>% 
  as.data.frame

genus.bin.dwr <- results(genus.bin.dwr) %>% 
  as.data.frame

genus.bin.som <- results(genus.bin.som) %>% 
  as.data.frame


##### HEATMAP NEGATIVE BINOMIAL MODELLING #####

library("pheatmap")

# Prepare data
genus.bin.ph.map  <- cbind(genus.norm, genus.bin.ph,  genus.hier)
genus.bin.den.map <- cbind(genus.norm, genus.bin.den, genus.hier)
genus.bin.dwr.map <- cbind(genus.norm, genus.bin.dwr, genus.hier)
genus.bin.som.map <- cbind(genus.norm, genus.bin.som, genus.hier)

# Keep only p < 0.05 and reorder by log2FoldChange
genus.bin.ph.map <- genus.bin.ph.map %>% 
  subset(pvalue < 0.1) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change)

genus.bin.den.map <- genus.bin.den.map %>% 
  subset(pvalue < 0.1) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change)

genus.bin.dwr.map <- genus.bin.dwr.map %>% 
  subset(pvalue < 0.1) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change)

genus.bin.som.map <- genus.bin.som.map %>% 
  subset(pvalue < 0.1) %>%
  mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
  arrange(Change)

# Remove unclassified
genus.bin.ph.map <- genus.bin.ph.map %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus))

genus.bin.den.map <- genus.bin.den.map %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus))

genus.bin.dwr.map <- genus.bin.dwr.map %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus))

genus.bin.som.map <- genus.bin.som.map %>% 
  subset(Genus != "") %>%
  subset(!grepl("Unclassified", Genus)) %>%
  subset(!grepl("Incertae Sedis", Genus))

# Plot
plotHeatmapBin(genus.bin.ph.map, "pH",
               filename = "02_METAXA/METAXA-pH-heatmap.png",
               annotation_row = genus.bin.ph.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.ph.map$Phylum,
                                               genus.bin.ph.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.den.map, "Density",
               filename = "02_METAXA/METAXA-den-heatmap.png",
               annotation_row = genus.bin.den.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.den.map$Phylum,
                                               genus.bin.den.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.dwr.map, "DryWetRatio",
               filename = "02_METAXA/METAXA-dwr-heatmap.png",
               annotation_row = genus.bin.dwr.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.dwr.map$Phylum,
                                               genus.bin.dwr.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.som.map, "SOM",
               filename = "02_METAXA/METAXA-SOM-heatmap.png",
               annotation_row = genus.bin.som.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.som.map$Phylum,
                                               genus.bin.som.map$Genus, sep = " / ")))

dev.off()

