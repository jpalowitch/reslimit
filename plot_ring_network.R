
library(ggplot2)
library(ggraph)
source("make_ring_network.R")
set.seed(12345)
layout1 <- create_layout(G1, layout = 'drl')
ggsave("ring.png",
       ggraph(layout1) + geom_edge_link() + geom_node_point())

layout2 <- create_layout(G2, layout = 'drl')
layout2$x <- layout1$x
layout2$y <- layout1$y

ggsave("mesh.png",
       ggraph(layout2) + geom_edge_link() + geom_node_point())