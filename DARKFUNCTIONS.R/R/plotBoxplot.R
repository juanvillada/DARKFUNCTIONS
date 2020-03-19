#' @export

plotBoxplot <- function(DATA, MAPPING) {
  ggplot(DATA, aes(x = !! enquo(MAPPING), y = Value)) +
    geom_boxplot(aes(fill = forcats::fct_rev(Layer)), outlier.shape = NA, position = position_dodge(preserve = "single")) +
    geom_point(aes(fill = forcats::fct_rev(Layer)), shape = 16, position = position_jitterdodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title = element_blank(),
          axis.line = element_blank(), panel.border = element_rect(fill = NA), strip.background = element_rect(fill = "grey90")) +
    scale_fill_manual(name = "Layer", values = c("#98c699", "#b7d8ff"))
}
