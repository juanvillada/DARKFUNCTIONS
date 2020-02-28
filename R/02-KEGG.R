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

# Transform to relative abundance
modules.rel <- modules %>% 
  mutate_all(function (x) x /sum(x))


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Compute richness
modules.rich <- computeRich(modules)

# Merge with metadata
modules.rich <- modules.rich %>% 
  full_join(metadata %>% select(Sample, Vegetation, Layer), by = "Sample")

# Plot
png("02_KEGG/KEGG-rich-boxplot.png", width = 1600, height = 1200, res = 150)
ggplot(modules.rich, aes(x = Vegetation, y = Richness)) +
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
modules.map <- bind_cols(modules.hier, modules.rel)

# Keep only pathways with relative abundance > 0.1%
modules.map <- modules.map %>% 
  filter(apply(modules.map[SAMPLES], 1, mean) > 0.001)

# Reorder taxa
modules.map <- modules.map %>% 
  arrange(desc(apply(modules.map[SAMPLES], 1, mean))) %>% 
  arrange(`Level 2`, `Level 3`)

# Transform to data.frame
modules.map <- data.frame(modules.map, row.names = rownames(modules.map), check.names = F)

# Plot
plotHeatmap(modules.map,
            filename = "02_KEGG/KEGG-heatmap.png",
            gaps_row = table(modules.map$`Level 2`) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Vegetation) %>% as.vector %>% cumsum,
            annotation_row = modules.map["Level 2"],
            annotation_colors = list(Vegetation = c(barren = "#dfc3f8", `deciduous shrub` = "#beefc1",
                                                    `evergreen shrub` = "#eca6c1", graminoid = "#61b7d9",
                                                    wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                   `Cellular processes` = "#feb8a6",
                                                   `Energy metabolism` = "#82e5ed",
                                                   `Environmental information processing` = "#e4b4e2",
                                                   `Genetic information processing` = "#a0eacd",
                                                   `Metabolism` = "#ffe0e7",
                                                   `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                   `Secondary metabolism` = "#aeac9b")),
            labels_row = as.character(paste(modules.map$`Level 3`,
                                            modules.map$`Level 4`, sep = " / ")))


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
modules.pmn <- adonis(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., 
                      filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), permutations = 9999, method = "bray")

# NMDS
modules.ord <- modules.rel %>% 
  t %>% 
  sqrt %>% 
  metaMDS(distance = "bray")

modules.ord.org <- modules.rel %>% 
  select(filter(metadata, Layer == "organic") %>% pull(Sample)) %>% 
  t %>% 
  sqrt %>% 
  metaMDS(distance = "bray")

# dbRDA
modules.rda <- ordiR2step(dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), distance = "bray"), 
                          dbrda(modules.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., filter(metadata, Sample %in% SAMPLES.FULL) %>% select(ORD.VARS), distance = "bray"), "forward")

# Plot
png("02_KEGG/KEGG-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(modules.ord, metadata)
while (!is.null(dev.list())) {
  dev.off()
}

png("02_KEGG/KEGG-NMDS-ORGANIC.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(modules.ord.org, metadata %>% filter(Layer == "organic"), 
                  ORDIELLIPSE = FALSE)
metadata %>% 
  filter(Layer == "organic") %>% 
  select(FLUX.VARS) %>% 
  envfit(modules.ord.org, .) %>% 
  plot(col = "#b3a5cb")
while (!is.null(dev.list())) {
  dev.off()
}

png("02_KEGG/KEGG-dbRDA_habitat.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(modules.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
while (!is.null(dev.list())) {
  dev.off()
}

png("02_KEGG/KEGG-dbRDA_layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(modules.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(modules.rda, display = "bp", col = "#e0a8d0")
while (!is.null(dev.list())) {
  dev.off()
}


##### NEGATIVE BINOMIAL MODELLING #####

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


##### HEATMAP NEGATIVE BINOMIAL MODELLING #####

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


##### FLUX DATA #####

SAMPLES.FULL2 <- metadata %>% 
  select(Sample, FLUX.VARS) %>% 
  drop_na %>% 
  pull(Sample) %>% 
  as.vector

# Run models
modules.bin.ch <- runDeseq(modules %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "CH4")
modules.bin.co <- runDeseq(modules %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "CO2")
modules.bin.no <- runDeseq(modules %>% select(SAMPLES.FULL2), metadata %>% filter(Sample %in% SAMPLES.FULL2), "N2O")

# Get normalised counts
modules.norm2 <- counts(modules.bin.ch, normalized = T) %>% 
  as_tibble

# Get results
modules.bin.ch <- results(modules.bin.ch) %>% 
  as_tibble

modules.bin.co <- results(modules.bin.co) %>% 
  as_tibble

modules.bin.no <- results(modules.bin.no) %>% 
  as_tibble

# Prepare data
modules.bin.ch.map <- bind_cols(modules.norm2, modules.bin.ch, modules.hier)
modules.bin.co.map <- bind_cols(modules.norm2, modules.bin.co, modules.hier)
modules.bin.no.map <- bind_cols(modules.norm2, modules.bin.no, modules.hier)

# Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
modules.bin.ch.map <- filterDeseq2(modules.bin.ch.map)
modules.bin.co.map <- filterDeseq2(modules.bin.co.map)
modules.bin.no.map <- filterDeseq2(modules.bin.no.map)

# Transform to data.frame
modules.bin.ch.map <- data.frame(modules.bin.ch.map, row.names = rownames(modules.bin.ch.map), check.names = F)
modules.bin.co.map <- data.frame(modules.bin.co.map, row.names = rownames(modules.bin.co.map), check.names = F)
modules.bin.no.map <- data.frame(modules.bin.no.map, row.names = rownames(modules.bin.no.map), check.names = F)

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

plotHeatmapBin2(modules.bin.ch.map, "CH4",
                filename = "02_KEGG/KEGG-CH4-heatmap.png",
                annotation_row = cbind(modules.bin.ch.map["Level 2"], modules.bin.ch.map["Change"]),
                annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                         `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Cellular processes` = "#feb8a6",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Environmental information processing` = "#e4b4e2",
                                                       `Genetic information processing` = "#a0eacd",
                                                       `Metabolism` = "#ffe0e7",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                       `Secondary metabolism` = "#aeac9b")),
                labels_row = as.character(paste(modules.bin.ch.map$`Level 3`,
                                                modules.bin.ch.map$`Level 4`, sep = " / ")))

plotHeatmapBin2(modules.bin.co.map, "CO2",
                filename = "02_KEGG/KEGG-CO2-heatmap.png",
                annotation_row = cbind(modules.bin.co.map["Level 2"], modules.bin.co.map["Change"]),
                annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                         `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Cellular processes` = "#feb8a6",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Environmental information processing` = "#e4b4e2",
                                                       `Genetic information processing` = "#a0eacd",
                                                       `Metabolism` = "#ffe0e7",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                       `Secondary metabolism` = "#aeac9b")),
                labels_row = as.character(paste(modules.bin.co.map$`Level 3`,
                                                modules.bin.co.map$`Level 4`, sep = " / ")))

plotHeatmapBin2(modules.bin.no.map, "N2O",
                filename = "02_KEGG/KEGG-N2O-heatmap.png",
                annotation_row = cbind(modules.bin.no.map["Level 2"], modules.bin.no.map["Change"]),
                annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd"),
                                         `Level 2` = c(`Carbohydrate and lipid metabolism` = "#7ec5ef",
                                                       `Cellular processes` = "#feb8a6",
                                                       `Energy metabolism` = "#82e5ed",
                                                       `Environmental information processing` = "#e4b4e2",
                                                       `Genetic information processing` = "#a0eacd",
                                                       `Metabolism` = "#ffe0e7",
                                                       `Nucleotide and amino acid metabolism` = "#94b2a6",
                                                       `Secondary metabolism` = "#aeac9b")),
                labels_row = as.character(paste(modules.bin.no.map$`Level 3`,
                                                modules.bin.no.map$`Level 4`, sep = " / ")))


