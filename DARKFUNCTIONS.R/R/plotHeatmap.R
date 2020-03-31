#' @export

plotHeatmap <- function(DATA, annotation_colors, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem, Layer) %>%
    pull(Sample) %>%
    as.vector

  ANNOTATION_COL <- data.frame(Layer = metadata %>% select(Layer),
                               Ecosystem = metadata %>% select(Ecosystem),
                               row.names = metadata %>% pull(Sample))

  ANNOTATION_COLORS <- c(list(Ecosystem = c(barren = "#dfc3f8", heathland = "#61b7d9", wetland = "#f9b99f"),
                              Layer = c(mineral = "#b7d8ff", organic = "#98c699")),
                         annotation_colors)

  DATA <- data.frame(DATA %>%
                       select(all_of(ORDER)) %>%
                       sqrt, row.names = rownames(DATA))

  pheatmap(DATA, color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
           annotation_col = ANNOTATION_COL, annotation_colors = ANNOTATION_COLORS, ...)
}
