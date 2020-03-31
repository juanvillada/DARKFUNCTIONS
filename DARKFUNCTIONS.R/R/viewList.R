#' @export

viewList <- function(x, y = NULL) {
  if (! is.null(y)) {
    cat(y, "\n")
    cat(rep("-", nchar(y)), "\n", sep = "")
    print(x[[y]])
  }
  else {
    for (y in names(x)) {
      cat(y, "\n")
      cat(rep("-", nchar(y)), "\n", sep = "")
      print(x[[y]])
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}
