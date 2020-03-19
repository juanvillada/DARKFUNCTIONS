#' @export

filterDeseq <- function(x, y) {
  x %>%
    filter(pvalue < y) %>%
    mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>%
    arrange(Change)
}
