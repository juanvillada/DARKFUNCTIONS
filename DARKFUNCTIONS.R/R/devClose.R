#' @export

devClose <- function() {
  while (!is.null(dev.list())) {
    dev.off()
  }
}
