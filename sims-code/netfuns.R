library(igraph)

twocliq <- function (nS, nC, p) {
  
  # Making graphs
  gS <- erdos.renyi.game(nS, p)
  gC1 <- gC2 <- erdos.renyi.game(nC, 1)
  
  # Getting disjoint edgelists
  eL <- get.edgelist(gS)
  eL1 <- get.edgelist(gC1) + nS
  eL2 <- get.edgelist(gC2) + nS + nC
  
  # Finding connector nodes
  minS <- min(eL); maxS <- max(eL)
  minC1 <- min(eL1); maxC1 <- max(eL1)
  minC2 <- min(eL2); maxC2 <- max(eL2)
  
  # Combining edgelists and adding connections
  eL <- rbind(eL, eL1, eL2, matrix(c(minS, minC1,
                                     maxS, minC2,
                                     maxC1, maxC2),
                                   ncol = 2,
                                   byrow = TRUE))
  
  G <- graph.edgelist(eL, directed = FALSE)
  return(G)
}

