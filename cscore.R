cscore <- function (G, C) {
  
  d <- degree(G)
  N <- length(V(G))
  m <- 2 * nrow(edges)
  
  # Setup for C
  dC <- d[C]
  edges <- get.edgelist(G)
  adj <- get.adjacency(G)
  kints <- colSums(adj[C, C])
  mC_int <- sum(kints)
  mC <- sum(dC)
  mC_ext <- mC - mC_int
  mstar <- m - mC
 
  # f-scores for C
  f_calc <- function (i) {
    dk <- dC[i]
    logf <- phyper(k, mC_ext, mstar, dk, log = TRUE, lower.tail = FALSE)
    return(sum(exp(logf)))
  }
  fscores <- phyper(kints - 1, mC_ext, mstar - mC_ext, dC, lower.tail = FALSE)
  worst_node <- C[which.max(fscores)]
  worst_node2 <- C[which.max(fscores[-worst_node])]
  
  # setup for C \ {w}
  C2 <- setdiff(C, worst_node)
  dC2 <- d[C2]
  kints2 <- kints[!C == worst_node]
  mC2_int <- sum(kints2)
  mC2 <- sum(dC2)
  mC2_ext <- mC2 - mC2_int
  mtild <- (N - length(C2)) * kints2[C2 == worst_node2]
  
  g_calc <- function (i) {
    dk <- dC2[i]
    logf <- sapply(kints2[i]:dk, function (k) dhyper(k, mC_ext, mstar, dk, log = TRUE))
    return(sum(exp(logf)))
  }
  
  
  
  
  
  
}