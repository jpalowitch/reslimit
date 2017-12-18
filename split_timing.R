if (FALSE) {
  edge.list.sub <- edge.list[edge.list[, 1] %in% comm, ]
  in.counts1 <- table(edge.list.sub[, 2])
  in.counts2 <- tapply(edge.list.sub[, 2], edge.list.sub[, 2], length)
  in.counts3 <- unlist(lapply(split(edge.list.sub[, 2], edge.list.sub[, 2]), length))
  
  
  in.counts2 <- colSums(adj.net[comm, ])
  make.adj <- graph.edgelist(edge.list, directed=FALSE)
}

nsims <- 30
ns <- 10^(2:7)
ps <- seq(0.1, 0.5, 0.1)

times <- array(dim=c(length(ns), length(ps), 3))
library(rbenchmark)

for (i in seq_along(ns)) {
  for (j in seq_along(ps)) {
  
    n <- ns[i]
    p <- ps[j]
    cat("n =", n, "p =", p, "\n")
    vec <- sample(n * p, n, replace=TRUE)
    
    times_ij <- benchmark(meth1=table(vec),
                          meth2=tapply(vec, vec, length),
                          meth3=unlist(lapply(split(vec, vec), length)),
                          replications=nsims)
    times[i, j, ] <- times_ij$elapsed
    
  }
}

library(reshape2)
library(ggplot2)
plotdf <- melt(times)
names(plotdf) <- c("n", "p", "Method", "Time")
plotdf$n <- ns[plotdf$n]
plotdf$p <- ps[plotdf$p]
plotdf$Time <- 1000 * plotdf$Time
plotdf$Method <- factor(plotdf$Method)
p <- ggplot(plotdf, aes(x=n, y=Time, colour=Method)) + 
  geom_line() + geom_point() + facet_wrap(~p) + 
  ggtitle("Vector Count Timing") + ylab("Milliseconds") + 
  scale_x_log10() + scale_y_log10()
p

