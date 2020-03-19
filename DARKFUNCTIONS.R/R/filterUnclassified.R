#' @export

filterUnclassified <- function(x, y) {
  x %>%
    filter(!!rlang::sym(y) != "") %>%
    filter(!grepl("Unclassified", !!rlang::sym(y))) %>%
    filter(!grepl("Incertae Sedis", !!rlang::sym(y)))
}
