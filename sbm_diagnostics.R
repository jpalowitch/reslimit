SBM_diagnostics <- function (G, truth_fn) {
  
  comms <- readLines(truth_fn)
  comms <- lapply(comms, function (C) as.numeric(strsplit(C, " ", fixed = TRUE)[[1]]))
  adjMat <- get.adjacency(G)
  K <- length(comms)
  avg_pout <- numeric(K)
  degs <- degree(G)
  
  for (i in 1:K) {
  
    nodes <- comms[[i]]
    pouts <- rowSums(adjMat[nodes, nodes]) / degs[nodes]
    avg_pout[i] <- mean(pouts)
    
  }
  
  return(list(avg_pout = avg_pout))
  
}