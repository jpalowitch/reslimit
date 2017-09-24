saveit <- function(..., file) {
    x <- list(...)
    save(list=names(x), file=file, envir=list2env(x))
}