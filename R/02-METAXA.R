library("tidyverse")
source("~/Scripts/DARKFUNCTIONS/R/00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Set working directory
setwd("~/Helsinki/analyses/")

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

# Read metadata
metadata <- readMetadata(SAMPLES)

# Read data
phylum <- read_delim("02_METAXA/summary_level_2.txt", delim = "\t")
genus <- read_delim("02_METAXA/summary_level_6.txt", delim = "\t")

# Split taxonomy
phylum <- phylum %>% 
  separate(Taxa, c("Domain", "Phylum"), ";")

genus <- genus %>% 
  separate(Taxa, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), ";")

# Keep only bacteria and archaea
phylum <- phylum %>% 
  filter(Domain %in% c("Bacteria", "Archaea"))

genus <- genus %>% 
  filter(Domain %in% c("Bacteria", "Archaea"))

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
  pull(Sample) %>% 
  as.vector

phylum <- phylum %>% 
  select(SAMPLES)

genus <- genus %>% 
  select(SAMPLES)

# Transform to relative abundance
phylum.rel <- phylum %>% 
  mutate_all(function (x) x /sum(x))

genus.rel <- genus %>% 
  mutate_all(function (x) x /sum(x))

# Summarise means and standard deviations
phylum.sum <- summariseMeans(phylum.rel, phylum.hier)
genus.sum <- summariseMeans(genus.rel, genus.hier)


##### BARPLOT #####

library("ggplot2")
library("forcats")

# Prepare data
phylum.bar <- bind_cols(Phylum = phylum.hier %>% select(Phylum), phylum.rel)

# Reorder taxa by abundance
phylum.bar <- phylum.bar %>% 
  arrange(desc(apply(phylum.bar %>% select(SAMPLES), 1, mean)))

# Melt and merge with metadata
phylum.bar <- phylum.bar %>% 
  gather(key = "Sample", value = "Count", -Phylum) %>% 
  full_join(metadata %>% select(Sample, Vegetation), by = "Sample") %>% 
  dplyr::rename(Taxa = Phylum) %>% 
  mutate(Sample = as.factor(Sample))

# Plot
png("02_METAXA/METAXA-PHYLUM-barplot.png", width = 1600, height = 1200, res = 150)
plot(plotBarplot(phylum.bar)) +
  facet_grid(cols = vars(Vegetation), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#ffbb9b", "#47e0ff", "#fbc387", "#64cbff", "#f0c683", "#6fb1ec", "#cfdd90", "#ce9cde", "#b9f5b5", "#faa0ca",
                               "#7ff9e1", "#ffaba1", "#59f0fc", "#d1a377", "#3ab9dc", "#ffddae", "#46d0c3", "#ffbdce", "#74dcb7", "#e9ccff",
                               "#91ba76", "#c3a1bf", "#e2ffb9", "#97acd4", "#73bc87", "#f0e3ff", "#68bb9e", "#d89d9a", "#c1ffdd", "#bcaa77",
                               "#bae8ff", "#fff4c2", "#75b8ae", "#eeffe2", "#91b3a0", "#c5fff0", "#a7af8c", "#b4e1e2"))
while (!is.null(dev.list())) {
  dev.off()
}


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
phylum.rich <- computeRich(phylum)
genus.rich <- computeRich(genus)

# Merge with metadata
genus.rich <- genus.rich %>% 
  full_join(metadata %>% select(Sample, Vegetation, Layer), by = "Sample")

# Plot
png("02_METAXA/METAXA-GENUS-rich-boxplot.png", width = 1600, height = 1200, res = 150)
ggplot(genus.rich, aes(x = Vegetation, y = Richness)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Vegetation)) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  facet_grid(rows = "Layer") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(),
        axis.line = element_blank(), legend.position = "none",
        panel.border = element_rect(fill = NA), strip.background = element_rect(fill = "grey90")) +
  ylab("Richness") +
  scale_fill_manual(values = c("#dfc3f8", "#beefc1", "#eca6c1", "#61b7d9", "#f9b99f", "#d8deff", "#b3a5cb"))
while (!is.null(dev.list())) {
  dev.off()
}


##### HEATMAP #####

library("pheatmap")

# Prepare data
phylum.map <- bind_cols(phylum.hier, phylum.rel)
genus.map <- bind_cols(genus.hier, genus.rel)

# For the phylum level, remove unclassified
phylum.map <- phylum.map %>% 
  filter(!grepl("Unclassified", Phylum))

# For the genus level, remove unclassified and keep only top 50 genera
genus.map <- genus.map %>% 
  arrange(desc(apply(genus.map[SAMPLES], 1, mean))) %>% 
  filter(Genus != "") %>%
  filter(!grepl("Unclassified", Genus)) %>%
  filter(!grepl("Incertae Sedis", Genus)) %>%
  head(50)

# Reorder taxa
phylum.map <- phylum.map %>% 
  arrange(Domain)

genus.map  <- genus.map %>% 
  arrange(Phylum, Genus)

# Transform to data.frame
phylum.map <- data.frame(phylum.map, row.names = rownames(phylum.map), check.names = F)
genus.map <- data.frame(genus.map, row.names = rownames(genus.map), check.names = F)

# Plot
plotHeatmap(phylum.map,
            filename = "02_METAXA/METAXA-PHYLUM-heatmap.png",
            gaps_row = table(phylum.map$Domain) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Vegetation) %>% as.vector %>% cumsum,
            annotation_row = phylum.map["Domain"],
            annotation_colors = list(Vegetation = c(barren = "#dfc3f8", `deciduous shrub` = "#beefc1",
                                                    `evergreen shrub` = "#eca6c1", graminoid = "#61b7d9",
                                                    wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     Domain = c(Archaea = "#8befff", Bacteria = "#acaf79")),
            labels_row = as.character(phylum.map$Phylum))

plotHeatmap(genus.map,
            filename = "02_METAXA/METAXA-GENUS-heatmap.png",
            gaps_row = table(genus.map$Phylum) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Vegetation) %>% as.vector %>% cumsum,
            annotation_row = genus.map["Phylum"],
            annotation_colors = list(Vegetation = c(barren = "#dfc3f8", `deciduous shrub` = "#beefc1",
                                                    `evergreen shrub` = "#eca6c1", graminoid = "#61b7d9",
                                                    wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                Verrucomicrobia = "#93b48d")),
            labels_row = as.character(genus.map$Genus))


##### ORDINATION #####

library("vegan")

# Create list of variables for ordination
ORD.VARS <- c("SoilWater", "pH", "SOM")
FLUX.VARS <- c("CH4", "CO2", "N2O")

# Create list of samples with full metadata
SAMPLES.FULL <- metadata %>% 
  select(Sample, ORD.VARS) %>% 
  drop_na %>% 
  pull(Sample) %>% 
  as.vector

# NP-MANOVA
genus.pmn <- adonis(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., 
                    filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), permutations = 9999, method = "bray")

# NMDS
genus.ord <- genus.rel %>% 
  t %>% 
  sqrt %>% 
  metaMDS(distance = "bray")

genus.ord.org <- genus.rel %>% 
  select(filter(metadata, Layer == "organic") %>% pull(Sample)) %>% 
  t %>% 
  sqrt %>% 
  metaMDS(distance = "bray")

# dbRDA
genus.rda <- ordiR2step(dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), distance = "bray"), 
                        dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), distance = "bray"), "forward")

# Plot
png("02_METAXA/METAXA-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(genus.ord, metadata)
while (!is.null(dev.list())) {
  dev.off()
}

png("02_METAXA/METAXA-NMDS-ORGANIC.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(genus.ord.org, metadata %>% filter(Layer == "organic"), 
                  ORDIELLIPSE = FALSE)
metadata %>% 
  filter(Layer == "organic") %>% 
  select(FLUX.VARS) %>% 
  envfit(genus.ord.org, .) %>% 
  plot(col = "#b3a5cb")
while (!is.null(dev.list())) {
  dev.off()
}

png("02_METAXA/METAXA-dbRDA_habitat.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(genus.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
while (!is.null(dev.list())) {
  dev.off()
}

png("02_METAXA/METAXA-dbRDA_layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(genus.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
while (!is.null(dev.list())) {
  dev.off()
}


##### NEGATIVE BINOMIAL MODELLING #####

library("DESeq2")

# Run models
genus.bin.ph <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "pH")
genus.bin.wt <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SoilWater")
genus.bin.om <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SOM")

# Get normalised counts
genus.norm <- counts(genus.bin.ph, normalized = T) %>% 
    as_tibble

# Get results
genus.bin.ph <- results(genus.bin.ph) %>% 
  as_tibble

genus.bin.wt <- results(genus.bin.wt) %>% 
  as_tibble

genus.bin.om <- results(genus.bin.om) %>% 
  as_tibble


##### HEATMAP NEGATIVE BINOMIAL MODELLING #####

library("pheatmap")

# Prepare data
genus.bin.ph.map <- bind_cols(genus.norm, genus.bin.ph, genus.hier)
genus.bin.wt.map <- bind_cols(genus.norm, genus.bin.wt, genus.hier)
genus.bin.om.map <- bind_cols(genus.norm, genus.bin.om, genus.hier)

# Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
genus.bin.ph.map <- filterDeseq(genus.bin.ph.map)
genus.bin.wt.map <- filterDeseq(genus.bin.wt.map)
genus.bin.om.map <- filterDeseq(genus.bin.om.map)

# Transform to data.frame
genus.bin.ph.map <- data.frame(genus.bin.ph.map, row.names = rownames(genus.bin.ph.map), check.names = F)
genus.bin.wt.map <- data.frame(genus.bin.wt.map, row.names = rownames(genus.bin.wt.map), check.names = F)
genus.bin.om.map <- data.frame(genus.bin.om.map, row.names = rownames(genus.bin.om.map), check.names = F)

# Plot
plotHeatmapBin(genus.bin.ph.map, "pH",
               filename = "02_METAXA/METAXA-pH-heatmap.png",
               annotation_row = genus.bin.ph.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.ph.map$Phylum,
                                               genus.bin.ph.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.wt.map, "SoilWater",
               filename = "02_METAXA/METAXA-SoilWater-heatmap.png",
               annotation_row = genus.bin.wt.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.wt.map$Phylum,
                                               genus.bin.wt.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.om.map, "SOM",
               filename = "02_METAXA/METAXA-SOM-heatmap.png",
               annotation_row = genus.bin.om.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.om.map$Phylum,
                                               genus.bin.om.map$Genus, sep = " / ")))


##### FLUX DATA #####

SAMPLES.FULL2 <- metadata %>% 
  select(Sample, FLUX.VARS) %>% 
  drop_na %>% 
  pull(Sample) %>% 
  as.vector

# Run models
genus.bin.ch <- runDeseq(genus %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "CH4")
genus.bin.co <- runDeseq(genus %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "CO2")
genus.bin.no <- runDeseq(genus %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "N2O")

# Get normalised counts
genus.norm2 <- counts(genus.bin.ch, normalized = T) %>% 
  as_tibble

# Get results
genus.bin.ch <- results(genus.bin.ch) %>% 
  as_tibble

genus.bin.co <- results(genus.bin.co) %>% 
  as_tibble

genus.bin.no <- results(genus.bin.no) %>% 
  as_tibble

# Prepare data
genus.bin.ch.map <- bind_cols(genus.norm2, genus.bin.ch, genus.hier)
genus.bin.co.map <- bind_cols(genus.norm2, genus.bin.co, genus.hier)
genus.bin.no.map <- bind_cols(genus.norm2, genus.bin.no, genus.hier)

# Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
genus.bin.ch.map <- filterDeseq(genus.bin.ch.map)
genus.bin.co.map <- filterDeseq(genus.bin.co.map)
genus.bin.no.map <- filterDeseq(genus.bin.no.map)

# Transform to data.frame
genus.bin.ch.map <- data.frame(genus.bin.ch.map, row.names = rownames(genus.bin.ch.map), check.names = F)
genus.bin.co.map <- data.frame(genus.bin.co.map, row.names = rownames(genus.bin.co.map), check.names = F)
genus.bin.no.map <- data.frame(genus.bin.no.map, row.names = rownames(genus.bin.no.map), check.names = F)

# Plot
plotHeatmapBin2 <- function (x, y, ...) {
  METADATA <- metadata %>% 
    filter(Sample %in% SAMPLES.FULL2)
  
  ORDER <- METADATA %>% 
    arrange(!! rlang::sym(y)) %>% 
    select(Sample) %>% 
    unlist %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(METADATA[y], row.names = METADATA$Sample)
  
  pheatmap(x[ORDER], color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           annotation_col = ANNOTATION_COL[y],
           gaps_row = subset(x, Change == "Negative") %>% nrow, ...)
}

plotHeatmapBin2(genus.bin.ch.map, "CH4",
               filename = "02_METAXA/METAXA-CH4-heatmap.png",
               annotation_row = genus.bin.ch.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.ch.map$Phylum,
                                               genus.bin.ch.map$Genus, sep = " / ")))

plotHeatmapBin2(genus.bin.co.map, "CO2",
                filename = "02_METAXA/METAXA-CO2-heatmap.png",
                annotation_row = genus.bin.co.map["Change"],
                annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
                labels_row = as.character(paste(genus.bin.co.map$Phylum,
                                                genus.bin.co.map$Genus, sep = " / ")))

plotHeatmapBin2(genus.bin.no.map, "N2O",
                filename = "02_METAXA/METAXA-N2O-heatmap.png",
                annotation_row = genus.bin.no.map["Change"],
                annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
                labels_row = as.character(paste(genus.bin.no.map$Phylum,
                                                genus.bin.no.map$Genus, sep = " / ")))

