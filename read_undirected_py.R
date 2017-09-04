library(igraph)

read_undirected_py <- function (fn) {
  el <- as.matrix(read.table(fn))
  N <- length(unique(as.vector(el)))
  order_seen <- as.vector(t(el))
  lookup <- integer(N)
  lookup[unique(order_seen)] <- 1:N
  order_seen2 <- lookup[order_seen]
  el2 <- matrix(order_seen2, ncol = 2, byrow = TRUE)
  swapme <- el2[ , 2] < el2[ , 1]
  el2[swapme, ] <- el2[swapme, 2:1]
  g2 <- graph.edgelist(el2, directed = FALSE)
  return(g2)
}