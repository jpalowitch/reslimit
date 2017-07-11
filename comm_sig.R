comm_sig <- function (membership, G) {
  
  
  K <- max(membership)
  d <- degree(G)
  m <- sum(d) / 2
  el <- get.edgelist(G)
  mem_el <- apply(el, 2, function (r) membership[r])
  pvals <- numeric(K)
  
  for (i in 1:K) {
    
    dCi <- 2 * sum(mem_el[ , 1] == i & mem_el[ , 2] == i)
    dC <- sum(d[membership == i])
    pvals[i] <- pbinom(dCi - 1, dC, dC / (2 * m), lower.tail = FALSE)
    
  }
  
  return(pvals)
  
}