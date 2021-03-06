source("full_dcsbm.R")
source("cluster_resolution.R")
source("recursive_mod.R")
source("igraph_recast.R")
source("plot_hierarchies.R")

N <- 1000

P <- matrix(c(0.8, 0.1, 0.2, 0.2,
              0.1, 0.8, 0.2, 0.2,
              0.2, 0.2, 0.6, 0.4,
              0.2, 0.2, 0.4, 0.6),
            ncol = 4)
membership <- rep(1:4, each = N / 4)
degrees <- rep(N / 2, N)
param_list <- make_param_list()
set.seed(12345)
test_net <- DCSBM(param_list, degrees = degrees, type = "slow",
                  P = P, membership = membership, adjust = FALSE)
edgelist1 <- get.edgelist(test_net$graph)
png("hierarchical_adj.png")
image(get.adjacency(test_net$graph))
dev.off()

test_net_el_recast <- igraph_recast(edgelist1, membership)

write.table(test_net_el_recast$edgelist, sep = "\t",
            file = "test_net.dat", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(as.character(test_net_el_recast$membership), 
           con = "test_net_membership.dat")

G <- graph.edgelist(as.matrix(read.table("test_net.dat")), directed = FALSE)
membership <- as.numeric(readLines("test_net_membership.dat"))
plotorder <- order(membership)
image(get.adjacency(G)[plotorder, plotorder])


res <- recursive_mod(fn = "test_net.dat", cscore_type = "r_bscore")

res_search <- cluster_resolution(fn = "test_net.dat", res_start = 0, res_end = 5, interval = 0.1,
                                 test_partitions = res$alllevels)
source("plot_hierarchies.R")
plot_hierarchies(res_search, fn = "test_net_res.png")

# Running SBM
write.graph(G, format = "graphml", file = "test_net.graphml")
curr_dir <- '.'
system(paste("/usr/bin/python",
             "fit_sbm.py",
             file.path(curr_dir, paste0("test_net", ".graphml")),
             file.path(curr_dir, paste0("test_net", "_gtMemship.dat"))))
sbm_membership <- as.integer(readLines(file.path(curr_dir, paste0("test_net", "_gtMemship.dat")))) + 1

# comparing rmod and truth
compare(res$alllevels[[2]], membership, method = "nmi")
compare(sbm_membership, membership, method = "nmi")
