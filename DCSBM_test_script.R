library(Matrix)
source("full_dcsbm.R")
param_list <- make_param_list(min_c = 20, max_c = 100, k = 20, max_k = 50, s2n = 10, N = 1000)
sbm_obj <- DCSBM(param_list, muversion = FALSE, type = "slow")
K <- max(sbm_obj$membership)
comms <- lapply(1:K, function (i) which(sbm_obj$membership == i))
nodeorder <- unlist(comms)
blockCounts <- matrix(0, K, K)
G <- sbm_obj$graph
#G <- graph.edgelist(as.matrix(read.table("sims-results/LFR_N=1000_B_20-100_increase_mu/1/1.dat")), directed = FALSE)
#comms <- readLines("sims-results/LFR_N=1000_B_20-100_increase_mu/1/1_truth.dat")
#comms <- lapply(comms, function (C) as.numeric(strsplit(C, split = " ")[[1]]))
#nodeorder <- unlist(comms)
adj <- get.adjacency(G)
K <- length(comms)
for (i in 1:K) {
  
  for (j in 1:K) {
    
    blockCounts[i, j] <- sum(adj[comms[[i]], comms[[j]]])
    
  }
  
}

good_nodes <- which(degree(G) > 0)
odegrees <- degree(G)[good_nodes]
pdegrees <- sbm_obj$degrees[good_nodes]
plot(pdegrees, odegrees)
abline(0, 1, col = "red")


image(get.adjacency(G)[nodeorder, nodeorder])
