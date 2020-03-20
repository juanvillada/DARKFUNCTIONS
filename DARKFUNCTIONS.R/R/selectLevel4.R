#' @export

selectLevel4 <- function(LEVEL4) {
  genes %>%
    right_join(.KO00002, by = c("level1", "level2", "level3", "level4", "KO")) %>%
    select(-gene) %>%
    inner_join(.KO00000, by = "KO") %>%
    filter(str_detect(level4, LEVEL4)) %>%
    mutate_at(MAGs, list(~replace(., is.na(.), 0))) %>%
    arrange(KO)
}
