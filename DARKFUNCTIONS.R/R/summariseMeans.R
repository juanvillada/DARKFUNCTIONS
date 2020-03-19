#' @export

summariseMeans <- function(DATA, METADATA) {
  bind_cols(METADATA,
            Mean = apply(DATA, 1, mean),
            SD = apply(DATA, 1, sd)) %>%
    arrange(desc(Mean)) %>%
    as_tibble
}
