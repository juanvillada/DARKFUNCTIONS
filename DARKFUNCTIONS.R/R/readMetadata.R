#' @export

readMetadata <- function(input) {
  input %>%
    read_delim(delim = "\t") %>%
    select(all_of(CAT.VARS), all_of(NUM.VARS), all_of(PROXY.VARS), all_of(FLUX.VARS)) %>%
    mutate(Sample = as.factor(Sample)) %>%
    mutate(Layer = as.factor(Layer)) %>%
    mutate(Ecosystem = as.factor(Ecosystem)) %>%
    mutate(Vegetation = as.factor(Vegetation)) %>%
    mutate_at(vars(all_of(NUM.VARS)), log) %>%
    mutate_at(vars(all_of(PROXY.VARS)), log)
}
