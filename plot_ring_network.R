
library(ggplot2)
library(ggraph)
source("make_ring_network.R")
set.seed(12345)
ring_networks <- make_ring_networks(10, nC = 6)
Gring <- ring_networks$ring; Gmesh <- ring_networks$mesh
#layout1 <- create_layout(Gring, layout = 'drl')
#ggsave("ring.png",
       #ggraph(layout1) + geom_edge_link() + geom_node_point())
#
#layout2 <- create_layout(Gmesh, layout = 'drl')
#layout2$x <- layout1$x
#layout2$y <- layout1$y
#
#ggsave("mesh.png",
       #ggraph(layout2) + geom_edge_link() + geom_node_point())

# Saving for graph-tool
library(igraph)
write.graph(Gring, file = "Gring.graphml", format = "graphml")
write.graph(Gmesh, file = "Gmesh.graphml", format = "graphml")
system("/usr/bin/python plot_rings.py")

