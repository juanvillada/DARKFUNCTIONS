#' @export

plotHeatmapKEGG2 <- function (DATA, RATIO) {
  p1 <- pheatmap(DATA$PRESENCE %>% t,
                 color = c("white", "orange"), silent = T,
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_row = paste(DATA$TAXONOMY %>% pull(Phylum), DATA$TAXONOMY %>% pull(MAG), sep = " | ") %>% as.character,
                 labels_col = paste(DATA$METADATA %>% pull(KO), DATA$METADATA %>% pull(gene), sep = " | ") %>% as.character)

  p2 <- pheatmap(DATA$COVERAGE %>% select(-Ecosystem) %>% sqrt %>% t,
                 color = colorRampPalette(colors = c("white", "yellow", "red"))(100), silent = T,
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_col = DATA$COVERAGE %>% pull(Ecosystem) %>% as.character,
                 show_rownames = F)

  plot_grid(as.grob(p2), as.grob(p1), ncol = 2, align = "h", rel_widths = RATIO)
}
