#' @export

extractGene <- function(MAG, KO) {
  KOtable <- lapply(kegg, function(x) {
    x[[MAG]] %>%
      getKOtable %>%
      filter(KO %in% !!KO) %>%
      select(-e_value) %>%
      unique
  }) %>%
    bind_rows %>%
    arrange(sequence) %>%
    unique

  sequences <- KOtable %>%
    pull(sequence) %>%
    unique

  gene_calls <- gene_calls %>%
    filter(gene_callers_id %in% sequences) %>%
    arrange(gene_callers_id)

  annotation <- gene_annot %>%
    filter(gene_callers_id %in% sequences) %>%
    filter(source != "COG_CATEGORY") %>%
    arrange(gene_callers_id)

  results <- list(KOtable = KOtable,
                  gene_calls = gene_calls,
                  annotation = annotation)

  return(results)
}
