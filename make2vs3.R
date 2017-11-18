source("full_dcsbm.R")
source("igraph_recast.R")
library(igraph)

N <- 300
set.seed(12345)
param_list <- make_param_list(N = N)
degrees <- rep(20, N)
membership <- c(rep(1, 2 * N / 3), rep(2, 1 * N / 6), rep(3, 1 * N / 6))
P <- matrix(0.1, 3, 3) + diag(3) * 0.7

sbm_obj <- DCSBM(param_list, degrees = degrees, membership = membership, P = P, type = "slow")
sbm_obj2 <- igraph_recast(as.matrix(get.edgelist(sbm_obj$graph)), membership = membership)
write.graph(graph.edgelist(sbm_obj2$edgelist, directed = FALSE),
            file = "2vs3.gml", format = "gml")
library(Matrix)
image(get.adjacency(sbm_obj$graph))
png("2vs3_adj.png")
image(get.adjacency(sbm_obj$graph))
dev.off()