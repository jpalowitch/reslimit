source("full_dcsbm.R")
source("cluster_resolution.R")
source("recursive_mod.R")
source("igraph_recast.R")

N <- 1000

P <- matrix(c(0.8, 0.1, 0.2, 0.2,
              0.1, 0.8, 0.2, 0.2,
              0.2, 0.2, 0.6, 0.4,
              0.2, 0.2, 0.4, 0.6),
            ncol = 4)
membership <- rep(1:4, each = N / 4)
degrees <- rep(N / 2, N)
param_list <- make_param_list()

test_net <- DCSBM(param_list, degrees = degrees, type = "slow",
                  P = P, membership = membership)
png("hierarchical_adj.png")
image(get.adjacency(test_net$graph))
dev.off()

test_net_el <- get.edgelist(test_net$graph)
test_net_el_recast <- igraph_recast(test_net_el)
membership <- membership[test_net_el_recast$lookup]
write.table(test_net_el_recast$edgelist, sep = "\t",
            file = "test_net.dat", quote = FALSE, row.names = FALSE, col.names = FALSE)


res <- recursive_mod(fn = "test_net.dat", cscore_type = "r_bscore")

res_search <- cluster_resolution(fn = "test_net.dat", res_start = 0, res_end = 5, interval = 0.1,
                                 test_partitions = res$alllevels)
source("plot_hierarchies.R")
plot_hierarchies(res_search, fn = "test_net_res.png")
