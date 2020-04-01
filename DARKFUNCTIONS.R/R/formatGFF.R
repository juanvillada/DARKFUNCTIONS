#' @export

formatGFF <- function(x, SOURCE) {
  gene_calls <- data.frame()
  for (y in names(x)) {
    gene_calls <- bind_rows(gene_calls,
                            x[[y]][["gene_calls"]])
  }

  annotation <- data.frame()
  for (y in names(x)) {
    annotation <- bind_rows(annotation,
                            x[[y]][["annotation"]])
  }

  annotation <- annotation %>%
    filter(source == SOURCE) %>%
    select(gene_callers_id, accession, `function`)

  gene_calls <- gene_calls %>%
    mutate(type = "CDS") %>%
    mutate(score = ".") %>%
    mutate(direction = ifelse(direction == "f", "+", "-")) %>%
    mutate(start = ifelse(start == 0, 1, start)) %>%
    mutate(phase = 0) %>%
    select(contig, source, type, start, stop, score, direction, phase, gene_callers_id)

  df <- left_join(gene_calls, annotation, by = "gene_callers_id") %>%
    separate(`function`, into = c("gene", "product"), sep = ";") %>%
    mutate(gene_callers_id = gsub("^", "ID=", gene_callers_id)) %>%
    mutate(gene = gsub("^", "Name=", gene)) %>%
    mutate(product = gsub("^", "product=", product)) %>%
    mutate(accession = gsub("^", "Alias=", accession)) %>%
    replace_na(list(accession = "", gene = "", product = "")) %>%
    unite(attributes, gene_callers_id, accession, gene, product, sep = ";")

  return(df)
}
