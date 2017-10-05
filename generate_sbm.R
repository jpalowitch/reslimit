suppressMessages(library(igraph))

generate_sbm <- function (membership, probs, out_degs = NULL,
                          python_path = "/usr/bin/python") {
  
  tmp <- tempdir()
  
  # Saving b
  write.table(membership - 1, quote = FALSE, 
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
  system(paste(python_path, "make_sbm.py", tmp))
  
  # Reading written graph
  G <- read.graph(file.path(tmp, "g.gml"), format = "gml")
  return(get.edgelist(G))

}