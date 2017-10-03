generate_sbm <- function (membership, probs, out_degs = NULL) {
  
  tmp <- tempdir()
  
  # Saving b
  write.table(membership, quote = FALSE, 
              row.names = FALSE, col.names = FALSE,
              file = file.path(tmp, "b.dat"))
  
  # Saving probs
  write.table(probs, quote = FALSE,
              row.names = FALSE, col.names = FALSE,
              file = file.path(tmp, "probs.dat"))
  
  # Saving out_degs
  if (!is.null(out_degs)) {
    write.table(out_degs, quote = FALSE, 
              row.names = FALSE, col.names = FALSE,
              file = file.path(tmp, "out_degs.dat"))
  }
  
  # Running python script
  
}