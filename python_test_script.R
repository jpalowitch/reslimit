library(igraph)
sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
sa_path <- "~/Documents/code/reslimit"
exper <- exper_names[2]
p <- 10
i <- 1

# Loading/setting parameters
rootdir <- file.path(sim_res_dir, exper)
load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
# Setting directory
curr_dir <- file.path(rootdir, p)
if (!dir.exists(curr_dir))
  dir.create(curr_dir)
# Naming the files
seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
dat_fn <- paste0(i, ".dat")

# Making R graph
#G <- graph.edgelist(as.matrix(read.table(file.path(curr_dir, dat_fn))), 
                          directed = FALSE)
fn <- file.path(curr_dir, dat_fn)

# Loading python stuff
require(rPython)
python.exec("import numpy")
python.exec("import louvain")
python.exec("import igraph as ig")

# Making python graph
python.exec(paste0("g = ig.Graph.Read_Ncol('", fn, "', directed = False)"))
python.exec("n = g.vcount()")
N <- python.get("n")

# Setting resolution parameter
res_par <- 1

# Running regular ouvain modularity
method <- 'Modularity'
python.exec(paste0("part = louvain.find_partition(g, method = '", method, "')"))
python.exec("membership = part.membership")
python.exec("mod = g.modularity(membership)")
mod_louvain <- python.get("mod")
mem_louvain <- python.get("membership")

# Running louvain modularity in R
source("read_undirected_py.R")
G <- read_undirected_py(fn)
mem_louvain_R <- cluster_louvain(G)$membership
modularity(G, mem_louvain_R)
modularity(G, mem_louvain + 1)

# Running RB modularity
python.exec(paste0("part = louvain.find_partition(g, method = '", "RBConfiguration", 
                   "', resolution_parameter = ", res_par, ")"))
python.exec("K = len(part)")
K <- python.get("K")
rbc_membership <- integer(N)
for (j in 1:K) {
  python.exec(paste0("x = part[", j - 1, "]"))
  rbc_membership[python.get("x")] <- j
}
if (sum(rbc_membership == 0) > 0) rbc_membership <- rbc_membership + 1

# Running cluster_louvain modularity in python
python.exec(paste0("part = louvain.find_partition(g, method = '", "RBConfiguration", 
                   "', resolution_parameter = ", res_par, ")"))
python.exec("K = len(part)")
K <- python.get("K")
rbc_membership <- integer(N)
for (j in 1:K) {
  python.exec(paste0("x = part[", j - 1, "]"))
  rbc_membership[python.get("x")] <- j
}
if (sum(rbc_membership == 0) > 0) rbc_membership <- rbc_membership + 1

head(mod_membership)
head(rbc_membership)
tail(mod_membership)
tail(rbc_membership)

mod_mod <- modularity(G, mod_membership)
rbc_mod <- modularity(G, rbc_membership)
