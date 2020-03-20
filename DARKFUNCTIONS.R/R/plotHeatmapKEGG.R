#' @export

plotHeatmapKEGG <- function (DATA, RATIO) {
  TAXONOMY <- taxonomy %>%
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species)

  MAGs <- TAXONOMY %>%
    pull(MAG)

  METADATA <- DATA %>%
    select(starts_with("level"))

  DATA <- DATA %>%
    select(all_of(MAGs))

  COVERAGE <- coverage.sum.eco %>%
    gather(key = "Ecosystem", value = value, 2:ncol(coverage.sum.eco)) %>%
    spread_(key = names(coverage.sum.eco)[1], value = "value") %>%
    select(Ecosystem, all_of(MAGs))

  p1 <- pheatmap(DATA %>% t,
                 color = c("white", "orange"), silent = T,
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_row = paste(TAXONOMY %>% pull(Phylum), TAXONOMY %>% pull(MAG), sep = " | ") %>% as.character,
                 labels_col = METADATA %>% pull(level4) %>% as.character)

  p2 <- pheatmap(COVERAGE %>% select(-Ecosystem) %>% sqrt %>% t,
                 color = colorRampPalette(colors = c("white", "yellow", "red"))(100), silent = T,
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_col = COVERAGE %>% pull(Ecosystem) %>% as.character, show_rownames = F)

  plot_grid(as.grob(p2), as.grob(p1), ncol = 2, align = "h", rel_widths = RATIO)
}
