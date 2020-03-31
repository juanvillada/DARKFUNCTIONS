#' @export

extractGene <- function(MAG, KO) {
  KOtable <- kegg[[MAG]] %>%
    getKOtable %>%
    filter(KO %in% !!KO) %>%
    arrange(sequence)

  sequences <- KOtable %>%
    pull(sequence) %>%
    unique

  gene_calls <- gene_calls %>%
    filter(gene_callers_id %in% sequences) %>%
    arrange(gene_callers_id)

  annotation <- gene_annot %>%
    filter(gene_callers_id %in% sequences) %>%
    arrange(gene_callers_id)

  results <- list(KOtable = KOtable,
                  gene_calls = gene_calls,
                  annotation = annotation)

  return(results)
}
