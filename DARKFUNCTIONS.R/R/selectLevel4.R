#' @export

selectLevel4 <- function (LEVEL4) {
  LEVEL4 <- KO00002 %>%
    filter(str_detect(`Level 4`, LEVEL4)) %>%
    pull(`Level 4`) %>%
    unique

  MAGs <- genes %>%
    filter(`Level 4` %in% LEVEL4) %>%
    select(all_of(MAGs)) %>%
    mutate_all(sum) %>%
    slice(1) %>%
    gather %>%
    filter(value > 0) %>%
    pull(key)

  TAXONOMY <- taxonomy %>%
    filter(MAG %in% MAGs) %>%
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species)

  MAGs <- TAXONOMY %>%
    pull(MAG)

  PRESENCE <- genes %>%
    filter(`Level 4` == LEVEL4) %>%
    select(all_of(MAGs))

  METADATA <- genes %>%
    filter(`Level 4` == LEVEL4) %>%
    select(starts_with("Level"), KO, gene)

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
