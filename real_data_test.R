library(igraph)
source("mod_sig.R")
source("../cscore_test/cscore/config_model.R")
source("recursive_mod.R")

### Bible
data_dir <- "data_sets"
data_set <- "bible"
data_fn  <- "names.dat"
edgelist <- read.table(file.path(data_dir, data_set, data_fn), skip = 5, sep = ',',
                       stringsAsFactors = FALSE)
all_names <- sort(unique(c(edgelist$V1, edgelist$V2)))
node1 <- match(edgelist$V1, all_names)
node2 <- match(edgelist$V2, all_names)
edgelist <- data.frame(node1, node2, edgelist$V3)
write.table(as.matrix(edgelist), quote = FALSE, row.names = FALSE, col.names = FALSE,
            file = file.path(data_dir, data_set, "edgelist.dat"))
res <- recursive_mod(file.path(data_dir, data_set, "edgelist.dat"), cscore_type = "r_bscore")


### Jazz
data_dir <- "data_sets"
data_set <- "jazz"
meta_fn  <- "meta.arenas-jazz"
data_fn  <- "out.arenas-jazz"
edgelist <- read.table(file.path(data_dir, data_set, data_fn), skip = 1, sep = "\t")
metadata <- readLines(file.path(data_dir, data_set, meta_fn))
sig_obj <- mod_sig(edgelist = edgelist)
write.table(as.matrix(edgelist), quote = FALSE, row.names = FALSE, col.names = FALSE,
            file = file.path(data_dir, data_set, "edgelist.dat"))
res <- recursive_mod(file.path(data_dir, data_set, "edgelist.dat"), cscore_type = "r_bscore")


### Copperfield
data_dir <- "data_sets"
data_set <- "copperfield"
data_fn  <- "adjnoun.gml"
copperfield <- read.graph(file.path(data_dir, data_set, data_fn), format = "gml")
edgelist <- get.edgelist(copperfield)
write.table(as.matrix(edgelist), quote = FALSE, row.names = FALSE, col.names = FALSE,
            file = file.path(data_dir, data_set, "edgelist.dat"))
res <- recursive_mod(file.path(data_dir, data_set, "edgelist.dat"), cscore_type = "r_bscore")
sig_obj <- mod_sig(edgelist = edgelist)

