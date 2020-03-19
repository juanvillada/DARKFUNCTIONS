#' @export

computeRich <- function(DATA) {
  tibble(Sample = SAMPLES,
         Value = apply(DATA, 2, function (y) sum(y > 0))) %>%
    full_join(metadata %>% select(all_of(CAT.VARS)), by = "Sample")
}
