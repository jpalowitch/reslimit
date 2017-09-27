make_ring_networks <- function (k, nC = 10, eB = 1) {
  
  suppressMessages(require(igraph))
  
  # Making cliques
  cdf <- t(combn(nC, 2))
  cdfs <- lapply(1:k, function (i) cdf + (i - 1) * nC)
  el <- do.call(rbind, cdfs)
  
  # Connecting some edges
  catch_edges <- nC * (0:(k - 1)) + 1
  pitch_edges <- nC * (1:k)
  pitch_edges <- c(pitch_edges[k], pitch_edges[1:(k - 1)])
  el1 <- rbind(el, cbind(pitch_edges, catch_edges))
  
  G1 <- graph.edgelist(el1, directed = FALSE)
  
  connections2 <- cbind(rep(pitch_edges, each = k),
                        rep(catch_edges, k))
  el2 <- rbind(el, connections2)
  G2 <- graph.edgelist(el2, directed = FALSE)
  return(list(ring = G1,
              mesh = G2))
  
}

