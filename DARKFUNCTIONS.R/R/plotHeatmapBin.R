#' @export

plotHeatmapBin <- function (DATA, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem) %>%
    pull(Sample) %>%
    as.vector

  ANNOTATION_COL <- data.frame(metadata %>% select(Ecosystem), row.names = metadata %>% pull(Sample))

  pheatmap(DATA %>% select(ORDER), color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
           annotation_col = ANNOTATION_COL,
           gaps_row = subset(DATA, Change == "Negative") %>% nrow, ...)
}
