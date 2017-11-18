# Set nsims
nsims <- 1000

# Set ns
Ns <- seq(150, 1200, 30)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

colPal <- gg_color_hue(2)

# set baseline prob and r
bp <- 0.2
r <- 0.7

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

mod2s <- mod3s <- matrix(0, length(Ns), nsims)

set.seed(12345)
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
    
    mod2s[j, i] <- mod2
    mod3s[j, i] <- mod3
    
  }
  
}

mod2means <- rowMeans(mod2s)
mod3means <- rowMeans(mod3s)
save(mod2means, mod3means, nsims, Ns, file = "sbm_mod_sim.RData")
#load("sbm_mod_sim.RData")
plotdf <- data.frame(n = rep(Ns, 2),
                     mods = c(mod2means, mod3means),
                     Partition = factor(rep(c("Joined", "Correct"), each = length(Ns))))
library(ggplot2)
my_theme <- theme(axis.text.x = element_text(size = 15),
                  axis.title.x = element_text(size = 15),
                  axis.text.y = element_text(size = 15),
                  axis.title.y = element_text(size = 15),
                  title = element_text(size = 18),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(size = 15))
p <- ggplot(plotdf, aes(x = n, y = mods, colour = Partition)) + geom_line() + 
  ylab("Average Modularity") + xlab("N") + 
  labs(title = "Correct Partition vs Joined Partition",
       subtitle = paste0(nsims, " repetitions per N")) + my_theme

ggsave("2vs3.png", p, width = 7, height = 7)
source("full_dcsbm.R")
source("igraph_recast.R")
library(igraph)

pd <- position_dodge(0.1) # move them .05 to the left and right

mod2s_t <- t(mod2s)
mod3s_t <- t(mod3s)
colnames(mod2s_t) <- colnames(mod3s_t) <- Ns
rownames(mod2s_t) <- rownames(mod3s_t) <- 1:nsims
library(reshape)
library(Rmisc)
plotdf2 <- melt(mod2s_t)
plotdf3 <- melt(mod3s_t)
plotdf2 <- cbind(plotdf2, part = "Joined")
plotdf3 <- cbind(plotdf3, part = "Correct")
plotdf_all <- rbind(plotdf2, plotdf3)
colnames(plotdf_all) <- c("Sample", "N", "Modularity", "Partition")
plotdf_all <- plotdf_all[-1]

tgc <- summarySE(plotdf_all, measurevar="Modularity", groupvars=c("N","Partition"))
tgc$N <- rep(Ns, each = 2)
pd <- position_dodge(0.1) # move them .05 to the left and right
p2 <- ggplot(tgc, aes(x=N, y=Modularity, colour=Partition)) + 
    geom_errorbar(aes(ymin=Modularity-ci, ymax=Modularity+ci), width=30, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) + 
  ylab("Average Modularity") + xlab("N") + 
  labs(title = "Correct Partition vs Joined Partition",
       subtitle = paste0(nsims, " repetitions per N")) + my_theme + 
  scale_colour_manual(values = rev(colPal))

ggsave("2vs3.png", p2, width = 7, height = 7)

N <- 150
set.seed(12345)
param_list <- make_param_list(N = N)
degrees <- rep(N / 2, N)
membership <- c(rep(1, 2 * N / 3), rep(2, 1 * N / 6), rep(3, 1 * N / 6))
P <- matrix(bp, 3, 3) + diag(3) * r

sbm_obj <- DCSBM(param_list, degrees = degrees, membership = membership, P = P, type = "slow")
sbm_obj2 <- igraph_recast(as.matrix(get.edgelist(sbm_obj$graph)), membership = membership)
write.graph(graph.edgelist(sbm_obj2$edgelist, directed = FALSE),
            file = "2vs3.gml", format = "gml")
library(Matrix)
image(get.adjacency(sbm_obj$graph))
png("2vs3_adj.png")
image(get.adjacency(sbm_obj$graph))
dev.off()