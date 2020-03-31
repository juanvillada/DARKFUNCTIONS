#' @export

plotBarplot <- function(DATA) {
  ggplot(DATA, aes(x = Sample, y = Count, fill = fct_reorder(Taxa, Count, mean))) +
    geom_bar(color = "black", size = 0.1, stat = "identity", position = "stack") +
    facet_grid(cols = vars(Ecosystem), scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),
          axis.line = element_blank(), legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(4, "pt"), legend.title = element_blank(), panel.border = element_rect(fill = NA),
          strip.background = element_rect(fill = "grey90")) +
    ylab("Percentage of mapped reads") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent)
}
