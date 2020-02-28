###  ### ### ### ### ### ### # ###
### DARKFUNCTIONS METAGENOMICS ###
### Auxiliary functions        ###
###  ### ### ### ### ### ### # ###


##### INITIALIZE ##### 

# Set working directory
setwd("~/Data/Helsinki/analyses/")

# Create list of variables
CAT.VARS <- c("Sample", "Layer", "Vegetation")
NUM.VARS <- c("SoilWater", "pH", "SOM")
PROXY.VARS <- c("SnowDepth", "Elevation")
FLUX.VARS <- c("CH4", "CO2", "N2O")


##### GENERAL ##### 

readMetadata <- function (x) {
  read_delim("00_METADATA/metadata.tsv", delim = "\t") %>% 
    select(Sample, Layer, Vegetation, SoilWater, pH, SOM, SnowDepth, Elevation, CH4, CO2, N2O) %>% 
    mutate(Sample = as.factor(Sample)) %>% 
    mutate(Layer = as.factor(Layer)) %>% 
    mutate(Vegetation = as.factor(Vegetation)) %>% 
    mutate_at(vars(NUM.VARS), log) %>% 
    mutate_at(vars(PROXY.VARS), log)
    filter(Sample %in% x)
}

runDeseq <- function (x, y, z) {
  x %>% 
    DESeqDataSetFromMatrix(., y, as.formula(paste("~", z))) %>%
    DESeq(fitType = "local", quiet = T)
}

filterDeseq <- function (x) {
  x %>% 
    filter(pvalue < 0.1) %>%
    mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
    arrange(Change) %>% 
    filter(Genus != "") %>%
    filter(!grepl("Unclassified", Genus)) %>%
    filter(!grepl("Incertae Sedis", Genus))
}

filterDeseq2 <- function (x) {
  x %>% 
    filter(pvalue < 0.05) %>%
    mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>% 
    arrange(Change, `Level 2`)
}



##### SUMMARIES ##### 

summariseMeans <- function (x, y) {
  data.frame(y,
             Mean = apply(x, 1, mean),
             SD = apply(x, 1, sd))
}

computeRich <- function (x) {
  tibble(Sample = names(x),
         Richness = apply(x, 2, function (y) sum(y > 0)))
}

##### PLOTTING #####

plotOrdinationLayer <- function (x, y, ORDIELLIPSE = TRUE, LEGEND = TRUE) {
  attach(y, warn.conflicts = F)
  
  COLORS <- c("#b7d8ff", "#98c699")
  
  plot(x, display = "sites", type = "n")
  points(x, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Layer])
  if (ORDIELLIPSE == TRUE) {
    ordiellipse(x, groups = Layer, kind = "sd", draw = "polygon", col = COLORS)
  }
  if (LEGEND == TRUE) {
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", col = COLORS, pch = c(16, 17))
  }
  detach(y)
}

plotOrdinationVeg <- function (x, y, ORDIELLIPSE = TRUE, LEGEND = TRUE) {
  attach(y, warn.conflicts = F)
  
  COLORS <- c("#dfc3f8", "#beefc1", "#eca6c1", "#61b7d9", "#f9b99f", "#d8deff")
  
  plot(x, display = "sites", type = "n")
  points(x, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Vegetation])
  if (ORDIELLIPSE == TRUE) {
    ordiellipse(x, groups = Vegetation, kind = "sd", draw = "polygon", col = COLORS)
  }
  if (LEGEND == TRUE) {
    legend("bottomright", legend = levels(metadata$Vegetation), bty = "n", fill = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
  }
  detach(y)
}

plotBarplot <- function (x) {
  ggplot(x, aes(x = Sample, y = Count, fill = fct_reorder(Taxa, Count, mean))) +
    geom_bar(color = "black", size = 0.1, stat = "identity", position = "stack") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),
          axis.line = element_blank(), legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(4, "pt"), legend.title = element_blank(), panel.border = element_rect(fill = NA),
          strip.background = element_rect(fill = "grey90")) +
    ylab("Percentage of mapped reads") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent)
}

plotHeatmap <- function (x, ...) {
  ORDER <- metadata %>% 
    arrange(Vegetation, Layer) %>% 
    pull(Sample) %>% 
    as.vector
  
  ANNOTATION_COL <- data.frame(Layer = metadata["Layer"], Vegetation = metadata["Vegetation"], row.names = metadata$Sample)
  
  pheatmap(sqrt(x[ORDER]), color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
           annotation_col = ANNOTATION_COL, ...)
  
  while (!is.null(dev.list())) {
    dev.off()
  }
}

plotHeatmapBin <- function (x, y, ...) {
  METADATA <- metadata %>% 
    filter(Sample %in% SAMPLES.FULL)
  
  ORDER <- METADATA %>% 
    arrange(!! rlang::sym(y)) %>% 
    select(Sample) %>% 
    unlist %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(METADATA[y], row.names = METADATA$Sample)
  
  pheatmap(x[ORDER], color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           annotation_col = exp(ANNOTATION_COL[y]),
           gaps_row = subset(x, Change == "Negative") %>% nrow, ...)
  
  while (!is.null(dev.list())) {
    dev.off()
  }
}
