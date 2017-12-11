# Getting software
library(ESSC)
library(reshape2)
system("bash get_CCME")
source("CCME/CCME.R")
sourceCpp("CCME/new_funs.cpp")

# Creating uniform 2D points
n <- 1000
set.seed(12345)
space_data <- matrix(runif(2 * n), ncol=2)

# Remove any points outside the circle
radii <- sqrt(rowSums((space_data - 0.5)^2))
space_data <- space_data[radii <= 0.5, , drop=FALSE]

# Finding distance matrix and .5-thresholded adjacency
space_data_dists <- space_data_adj <- as.matrix(dist(space_data))

xdists <- outer(space_data[ , 1], space_data[ , 1], "-")
ydists <- outer(space_data[ , 2], space_data[ , 2], "-")
pdists <- sqrt(xdists^2 + ydists^2)

space_data_adj[space_data_dists <= 0.5] <- 1
space_data_adj[space_data_dists > 0.5] <- 0

# Running ESSC
set.seed(12345)
essc_res <- essc(space_data_adj, alpha=0.05, verbose=TRUE, exhaustive=FALSE)

# Getting weighted edge_list
edge_list <- melt(space_data_dists)
names(edge_list) <- c("node1", "node2", "weight")
rownames(edge_list) <- NULL
edge_list <- edge_list[edge_list$node1 >= edge_list$node2, ]

# Weights are currently distances. Transforming to make into similarities
edge_list$weight <- max(edge_list$weight) - edge_list$weight

# Running CCME
set.seed(12345)
ccme_res <- CCME(edge_list)

# Plotting essc
K1 <- length(essc_res$Communities)
png("essc_comms.png", height=1000, width=ceiling(K2 / 2) * 500)
par(mfrow=c(2, ceiling(K2 / 2)))
for (i in 1:K1) {
  comm <- essc_res$Communities[[i]]
  plot(space_data[comm, ], col=i + 1, pch=i, cex=3, xlim=c(0, 1), ylim=c(0, 1),
       main=paste0("comm_", i), xlab="x", ylab="y")
  points(space_data[-comm, ], cex=2, pch=16)
}
dev.off()

# Plotting ccme
K2 <- length(ccme_res$communities)
png("ccme_comms.png", height=1000, width=ceiling(K2 / 2) * 500)
par(mfrow=c(2, ceiling(K2 / 2)))
for (i in 1:K2) {
  comm <- ccme_res$communities[[i]]
  plot(space_data[comm, ], col=i + 1, pch=i, cex=3, xlim=c(0, 1), ylim=c(0, 1),
       main=paste0("comm_", i), xlab="x", ylab="y")
  points(space_data[-comm, ], cex=2, pch=16)
}
dev.off()

