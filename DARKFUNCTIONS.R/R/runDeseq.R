#' @export

runDeseq <- function(x, y, z) {
  x %>%
    DESeqDataSetFromMatrix(., y, as.formula(paste("~", z))) %>%
    DESeq(fitType = "local", quiet = T)
}
