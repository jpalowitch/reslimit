source("recursive_mod.R")
source("cluster_resolution.R")
source("plot_hierarchies.R")
source("applications-code/plotPortComms2.R")

yearNum <- 2014
curr_dir <- file.path("applications-results", "airports", yearNum)
fn <- "year"
load(file.path(curr_dir, paste0(fn, ".RData")))
codeList <- rownames(data)

colnames(data) <- NULL
rownames(data) <- NULL
G <- graph.adjacency(data, mode = "undirected", weighted = "weight")

edge_list <- get.edgelist(G)
edge_list <- cbind(edge_list, E(G)$weight)
edge_list <- as.data.frame(edge_list)
names(edge_list) <- c("node1", "node2", "weight")
edge_list_unweighted <- edge_list[ , 1:2]

dat_fn <- paste0(yearNum, "_", fn, ".dat")


write.table(edge_list_unweighted,
            file = file.path(curr_dir, dat_fn),
            row.names = FALSE,
            col.names = FALSE)

res <- recursive_mod(file.path(curr_dir, dat_fn), cscore_type = "r_bscore")
