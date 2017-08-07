# Set path to compare
compare_path <- "cb-signi"

# Set network/comm path and names
net_path <- "sims2/er_increase_p/5"
net_fn <- "1"
com_fn <- "comm1"

# Calculating c-score with compare
netfile <- paste0(net_fn, ".dat")
comfile <- paste0(com_fn, ".dat")
system(paste(file.path(compare_path, 'compare'),
             '-f', file.path(net_path, netfile), 
             '-c', file.path(net_path, comfile),
             '-t 0.01 -cscore'))


# Reading c-score
cbscores <- read.table(file.path(net_path, paste0(net_fn, ".dat.table")), 
                       sep = "",
                       header = FALSE, row.names = 1)
cbscores <- unname(as.vector(cbscores))

# Calculating the c-score
edges <- as.matrix(read.table(netfile))
library(igraph)
library(Matrix)
G <- graph.edgelist(edges, directed = FALSE)
C <- as.numeric(strsplit(readLines(comfile), " ")[[1]])
