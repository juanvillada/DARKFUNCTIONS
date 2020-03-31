#' @export

selectLevel4 <- function(LEVEL4) {
  GENES <- genes %>%
    right_join(.KO00002, by = c("level1", "level2", "level3", "level4", "KO")) %>%
    select(-gene) %>%
    inner_join(.KO00000, by = "KO") %>%
    filter(str_detect(level4, LEVEL4)) %>%
    mutate_at(MAGs, list(~replace(., is.na(.), 0))) %>%
    arrange(KO)

  MAGs <- GENES %>%
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

  GENES_HIER <- GENES %>%
    select(starts_with("level"), KO, gene)

  GENES <- GENES %>%
    select(all_of(MAGs))

  ABUNDANCE <- bind_cols(taxonomy, abundance.rel) %>%
    select(MAG, all_of(SAMPLES)) %>%
    filter(MAG %in% MAGs) %>%
    arrange(match(MAG, MAGs))

  # ABUNDANCE.ECO <- abundance.eco %>%
  #   gather(key = "Ecosystem", value = value, 2:ncol(abundance.eco)) %>%
  #   spread_(key = names(abundance.eco)[1], value = "value") %>%
  #   select(Ecosystem, all_of(MAGs))

  results <- list(MAGs = MAGs,
                  TAXONOMY = TAXONOMY,
                  GENES = GENES,
                  GENES_HIER = GENES_HIER,
                  ABUNDANCE = ABUNDANCE)

  return(results)
}
