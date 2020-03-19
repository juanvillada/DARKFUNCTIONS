#' @export

selectLevel3 <- function (LEVEL3) {
  MODULES <- modules %>%
    filter(level3 == LEVEL3) %>%
    pull(level4)

  TAXONOMY <- taxonomy %>%
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species)

  MAGs <- TAXONOMY %>%
    pull(MAG)

  PRESENCE <- modules %>%
    filter(level4 %in% MODULES) %>%
    select(all_of(MAGs))

  METADATA <- modules %>%
    filter(level4 %in% MODULES) %>%
    select(starts_with("level"))

  COVERAGE <- coverage.sum.eco %>%
    gather(key = "Ecosystem", value = value, 2:ncol(coverage.sum.eco)) %>%
    spread_(key = names(coverage.sum.eco)[1], value = "value") %>%
    select(Ecosystem, all_of(MAGs))

  TAXONOMY <- data.frame(TAXONOMY, row.names = rownames(TAXONOMY), check.names = F)
  PRESENCE <- data.frame(PRESENCE, row.names = rownames(PRESENCE), check.names = F)
  COVERAGE <- data.frame(COVERAGE, row.names = rownames(COVERAGE), check.names = F)
  METADATA <- data.frame(METADATA, row.names = rownames(METADATA), check.names = F)

  return(list(TAXONOMY = TAXONOMY,
              PRESENCE = PRESENCE,
              METADATA = METADATA,
              COVERAGE = COVERAGE))
}
