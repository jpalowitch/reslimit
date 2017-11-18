# Set nsims
nsims <- 10

# Set ns
Ns <- seq(150, 1200, 30)

# set baseline prob and r
bp <- 0.2
r <- 0.6

# Set P
P1 <- matrix(c(0.06, 0.04, 0.00,
               0.04, 0.12, 0.04,
               0.00, 0.04, 0.66),
             ncol = 3)

P2 <- matrix(bp, ncol = 3, nrow = 3)
diag(P2) <- diag(P2) + r

tmp <- tempdir()
source("generate_sbm.R")
source("cluster_resolution.R")
source("igraph_recast.R")
source("full_dcsbm.R")

mod2s <- mod3s <- numeric(length(Ns))


for (j in seq_along(Ns)) {
  
  n <- Ns[j]
  
  for (i in 1:nsims) {
    
    cat("n =", n, "i =", i, "\n")
  
    # Making network
    commsizes <- n * c(2 / 3, 1 / 6, 1 / 6)
    membership <- rep.int(1:3, times = commsizes)
    probs <- P2 * n^2 / 9
    param_list <- make_param_list()
    sbm_sim <- DCSBM(param_list, type = "slow", degrees = rep(n / 2, n),
                     P = P2, membership = membership, adjust = FALSE)
    membership0 <- rep.int(1:2, times = c(commsizes[1], sum(commsizes[2:3])))
    
    # Calculating mods
    mod2 <- modularity(sbm_sim$graph, membership0)
    mod3 <- modularity(sbm_sim$graph, membership)
    
    mod2s[j] <- mod2s[j] + mod2 / nsims
    mod3s[j] <- mod3s[j] + mod3 / nsims
    
  }
  
}

plotdf <- data.frame(n = rep(Ns, 2),
                     mods = c(mod2s, mod3s),
                     partition = factor(rep(2:3, each = length(Ns))))
library(ggplot2)
p <- ggplot(plotdf, aes(x = n, y = mods, colour = partition)) + geom_line() + 
  ylab("Average Modularity, 10 repetitions") + xlab("N") + 
  ggtitle("Empirical Modularity of Correct (3) vs Joined (2) Partitions")

ggsave("2vs3.png", p, width = 7, height = 7)