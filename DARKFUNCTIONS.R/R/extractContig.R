#' @export

extractContig <- function(CONTIG) {
  gene_calls <- gene_calls %>%
    filter(contig %in% CONTIG) %>%
    arrange(gene_callers_id)

  sequences <- gene_calls %>%
    pull(gene_callers_id)

  annotation <- gene_annot %>%
    filter(gene_callers_id %in% sequences) %>%
    arrange(gene_callers_id)

  results <- list(gene_calls = gene_calls,
                  annotation = annotation)

  return(results)
}
