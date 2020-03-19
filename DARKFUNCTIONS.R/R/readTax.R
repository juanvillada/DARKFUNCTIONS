#' @export

readTax <- function() {
  bind_rows(read_delim("05_MAGs/GTDB/gtdbtk.bac120.summary.tsv", delim = "\t") %>%
              select(user_genome, classification) %>%
              rename(MAG = user_genome),
            read_delim("05_MAGs/GTDB/gtdbtk.ar122.summary.tsv", delim = "\t") %>%
              select(user_genome, classification) %>%
              rename(MAG = user_genome)) %>%
    filter(MAG %in% MAGs) %>%
    arrange(match(MAG, MAGs)) %>%
    mutate(classification = gsub("[a-z]__", "", classification)) %>%
    separate(classification, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
}
