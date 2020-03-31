#' @export

plotModules <- function(DATA, SCALING = 10, ...) {
  # p2 <- pheatmap(DATA$ABUNDANCE.ECO %>% select(-Ecosystem) %>% sqrt,
  #                color = colorRampPalette(colors = c("white", "yellow", "red"))(100), silent = T,
  #                cellwidth = SCALING, fontsize_row = SCALING, fontsize_col = SCALING,
  #                cluster_cols = F, cluster_rows = F, legend = F,
  #                labels_row = DATA$ABUNDANCE.ECO %>% pull(Ecosystem) %>% as.character, show_colnames = F, ...)
  #
  # plot(plot_grid(as.grob(p2), as.grob(p1), nrow = 2, align = "v", rel_heights = c(1, 10), greedy = T))

  pheatmap(DATA$GENES,
           color = c("white", "orange"),
           cellwidth = SCALING, cellheight = SCALING, fontsize_row = SCALING, fontsize_col = SCALING,
           cluster_cols = F, cluster_rows = F, legend = F,
           labels_row = paste(DATA$GENES_HIER %>% pull(KO), DATA$GENES_HIER %>% pull(gene), sep = " | ") %>% as.character,
           labels_col = paste(seq(1, length(DATA$MAGs)), DATA$TAXONOMY %>% pull(Phylum), DATA$TAXONOMY %>% pull(MAG), sep = ", ") %>% as.character)

  ORDER <- metadata %>%
    arrange(Ecosystem, Layer) %>%
    pull(Sample) %>%
    as.vector

  pheatmap(DATA$ABUNDANCE %>% select(all_of(ORDER)) %>% sqrt %>% t,
           color = colorRampPalette(colors = c("white", "yellow", "red"))(100),
           cellwidth = SCALING, cellheight = 7, fontsize_row = 7, fontsize_col = SCALING,
           cluster_cols = F, cluster_rows = F,
           annotation_row =  data.frame(Layer = metadata %>% select(Layer), Ecosystem = metadata %>% select(Ecosystem), row.names = metadata %>% pull(Sample)),
           annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#61b7d9", wetland = "#f9b99f"),
                                    Layer = c(mineral = "#b7d8ff", organic = "#98c699")),
           annotation_names_row = F, annotation_names_col = F,
           labels_col = seq(1, length(DATA$MAGs)) %>% as.character, ...)

  ABUNDANCE <- DATA$ABUNDANCE %>%
    gather(key= "Sample", value = "Abundance", -MAG) %>%
    group_by(Sample) %>%
    summarise(Value = sum(Abundance)) %>%
    full_join(metadata %>% select(all_of(CAT.VARS)))

  plotBoxplot(ABUNDANCE, Ecosystem)
}
