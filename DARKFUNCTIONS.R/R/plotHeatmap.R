#' @export

plotHeatmap <- function(DATA, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem, Layer) %>%
    pull(Sample) %>%
    as.vector

  ANNOTATION_COL <- data.frame(Layer = metadata %>% select(Layer), Ecosystem = metadata %>% select(Ecosystem), row.names = metadata %>% pull(Sample))

  pheatmap(DATA %>% select(ORDER) %>% sqrt, color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
           annotation_col = ANNOTATION_COL, ...)
}
