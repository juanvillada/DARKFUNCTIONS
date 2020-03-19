#' @export

plotHeatmapBin2 <- function (DATA, MAPPING, ...) {
  METADATA <- metadata %>%
    filter(Layer == "organic")

  ORDER <- METADATA %>%
    arrange(!! rlang::sym(MAPPING)) %>%
    pull(Sample) %>%
    as.vector

  ANNOTATION_COL <- data.frame(METADATA %>% select(MAPPING), row.names = METADATA %>% pull(Sample))

  pheatmap(DATA %>% select(ORDER), color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           annotation_col = ANNOTATION_COL,
           gaps_row = subset(DATA, Change == "Negative") %>% nrow, ...)
}
