###  ### ### ### ### ### ### # ###
### DARKFUNCTIONS METAGENOMICS ###
### Auxiliary functions        ###
###  ### ### ### ### ### ### # ###


##### TRANSFORMATIONS #####

metadataSub <- function (x) {
  subset(metadata, Sample %in% x)
}

metadataMerge <- function (x) {
  melt(x) %>% 
    merge(metadata %>% select(Sample, Habitat, Layer), by.x = "variable", by.y = "Sample")%>%
    setNames(c("Sample", "Taxa", "Count", "Habitat", "Layer"))
}


##### SUMMARIES ##### 

summariseMeans <- function (x, y) {
  data.frame(y,
             Mean = apply(x, 1, mean),
             SD = apply(x, 1, sd))
}

##### PLOTTING #####

plotOrdinationLayer <- function (x, y) {
  attach(y, warn.conflicts = F)
  
  COLORS <- c("#b7d8ff", "#98c699")
  
  plot(x, display = "sites", type = "n")
  points(x, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Layer])
  ordiellipse(x, groups = Layer, kind = "sd", draw = "polygon", col = COLORS)
  ordispider(x, groups = Layer, col = COLORS)
  legend("bottomleft", legend = levels(metadata$Layer), bty = "n", col = COLORS, pch = c(16, 17))
  
  detach(y)
}

plotOrdinationHabitat <- function (x, y) {
  attach(y, warn.conflicts = F)
  
  COLORS <- c("#dfc3f8", "#beefc1", "#eca6c1", "#61b7d9", "#f9b99f", "#d8deff","#b3a5cb")
  
  plot(x, display = "sites", type = "n")
  points(x, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Habitat])
  ordiellipse(x, groups = Habitat, kind = "sd", draw = "polygon", col = COLORS)
  ordispider(x, groups = Habitat, col = COLORS)
  legend("bottomright", legend = levels(metadata$Habitat), bty = "n", fill = COLORS)
  legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
  
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
  METADATA <- metadataSub(SAMPLES)
  
  ORDER <- METADATA[order(METADATA$Habitat, METADATA$Layer), ] %>% 
    select(Sample) %>% 
    unlist %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(Layer = METADATA["Layer"], Habitat = METADATA["Habitat"], row.names = METADATA$Sample)
  
  pheatmap(sqrt(x[ORDER]), color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
           annotation_col = ANNOTATION_COL, ...)
}

plotHeatmapBin <- function (x, y, ...) {
  METADATA <- metadataSub(SAMPLES.FULL)
  
  ORDER <- METADATA[order(METADATA[y]), ] %>% 
    select(Sample) %>% 
    unlist %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(METADATA[y], row.names = METADATA$Sample)
  
  pheatmap(x[ORDER], color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           annotation_col = exp(ANNOTATION_COL[y]),
           gaps_row = subset(x, Change == "Negative") %>% nrow, ...)
}
