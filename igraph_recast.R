library(igraph)

igraph_recast <- function (el, membership = NULL) {
  
  N <- length(unique(as.vector(el)))
  order_seen <- as.vector(t(el))
  lookup <- integer(N)
  lookup[unique(order_seen)] <- 1:N
  order_seen2 <- lookup[order_seen]
  el2 <- matrix(order_seen2, ncol = 2, byrow = TRUE)
  swapme <- el2[ , 2] < el2[ , 1]
  el2[swapme, ] <- el2[swapme, 2:1]
  
  # Getting new membership if it was put in
  if (!is.null(membership)) {
    membership_recast <- membership[order(lookup)]
  } else {
    membership_recast <- NULL
  }
  
  return(list(edgelist = el2, lookup = lookup, 
              membership = membership_recast))
  
}