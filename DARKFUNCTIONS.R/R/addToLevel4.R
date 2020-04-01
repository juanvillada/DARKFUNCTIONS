#' @export

addToLevel4 <- function(MODULE, KO_SELECT) {
  MAGs <- MODULE[["MAGs"]]

  GENES <- lapply(kegg, function(x) {
    lapply(MAGs, function(y) {
      df <- x[[y]] %>%
        getKOtable %>%
        select(KO) %>%
        filter(KO %in% KO_SELECT) %>%
        unique %>%
        inner_join(.KO00000, by = "KO") %>%
        mutate(MAG = y) %>%
        mutate(count = 1) %>%
        arrange(KO)

      if (nrow(df) == 0) {
        df <- tibble(KO = NA, gene = NA, MAG = y, count = 0)
      }

      return(df)
    }) %>%
      bind_rows %>%
      spread(MAG, count, fill = 0) %>%
      filter(!KO %in% NA) %>%
      select(KO, gene, all_of(MAGs))
  }) %>%
    bind_rows %>%
    group_by(KO, gene) %>%
    summarise_all(funs(sum(., na.rm = T))) %>%
    ungroup %>%
    mutate_at(MAGs, function (x) ifelse(x > 0, 1, 0))

  LEVEL1 <- MODULE[["GENES_HIER"]] %>% pull(level1) %>% unique
  LEVEL2 <- MODULE[["GENES_HIER"]] %>% pull(level2) %>% unique
  LEVEL3 <- MODULE[["GENES_HIER"]] %>% pull(level3) %>% unique
  LEVEL4 <- MODULE[["GENES_HIER"]] %>% pull(level4) %>% unique

  GENES <- GENES %>%
    mutate(level1 = LEVEL1, level2 = LEVEL2, level3 = LEVEL3, level4 = LEVEL4) %>%
    select(starts_with("level"), KO, gene, all_of(MAGs))

  GENES_OLD <- bind_cols(MODULE[["GENES_HIER"]], MODULE[["GENES"]])

  GENES <- bind_rows(GENES_OLD, GENES) %>%
    arrange(KO)

  GENES_HIER <- GENES %>%
    select(-all_of(MAGs))

  GENES <- GENES %>%
    select(all_of(MAGs))

  results <- list(MAGs = MAGs,
                  TAXONOMY = MODULE[["TAXONOMY"]],
                  GENES = GENES,
                  GENES_HIER = GENES_HIER,
                  ABUNDANCE = MODULE[["ABUNDANCE"]])

  return(results)
}
