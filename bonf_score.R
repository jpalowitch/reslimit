library(Matrix)

bonf_score <- function (group_s, edges, version = 1, borderp = 0.25, nsims = 100,
                        return_mean = FALSE) {
  
  if (length(group_s) < 3)
    return(1)
  
  # Finding degrees
  N <- max(edges)
  d1 <- unlist(lapply(split(edges[, 1], edges[, 1]), length))
  d2 <- unlist(lapply(split(edges[, 2], edges[, 2]), length))
  d <- integer(N)
  d1.names <- as.integer(names(d1))
  d2.names <- as.integer(names(d2))
  d[d1.names] <- d1
  d[d2.names] <- d[d2.names] + d2
  
  # Computing other basic values
  nodes <- 1:N
  m <- 2 * nrow(edges)
  edges.comm <- edges[edges[, 1] %in% group_s, ]
  in.counts <- unlist(lapply(split(edges.comm[, 2], edges.comm[, 2]), length))
  kints <- integer(N)
  kints[as.integer(names(in.counts))] <- in.counts
  kexts <- d - kints
  mC_int <- sum(kints[group_s])
  mC <- sum(d[group_s])
  mC_ext <- mC - mC_int
  mstar <- m - mC
  
  # Getting r-scores
  mC_ext_new <- mC_ext - kexts[group_s] + kints[group_s]
  mC_new <- mC - d[group_s]
  mstar_new <- m - mC_new
  r_scores <- phyper(kints[group_s] - 1, mC_ext_new, mstar_new - mC_ext_new,
                     d[group_s], lower.tail = FALSE)
  worst_node <- group_s[which.max(r_scores)]
  kw <- kints[worst_node]
  degree_ref <- d[worst_node]
    
  cs <- numeric(nsims)
  for (counter in 1:nsims) {
    
    # Getting random r-scores
    group_sc <- setdiff(nodes, group_s)
    r_scoresU <- phyper(kints[group_s] - 1, mC_ext_new, mstar_new - mC_ext_new,
                        d[group_s], lower.tail = FALSE)
    r_scoresL <- phyper(kints[group_s],     mC_ext_new, mstar_new - mC_ext_new,
                        d[group_s], lower.tail = FALSE)
    rand_rscores <- runif(length(group_s), r_scoresL, r_scoresU)
    
    # Doing c-score if version 1
    if (version == 1) {
      r1 <- sort(rand_rscores, decreasing = TRUE)[1]
      r2 <- sort(rand_rscores, decreasing = TRUE)[2]
      cs[counter] <- 1 - (1 - (r1 - r2) / (1 - r2))^(N - length(group_s) + 1)
    }
    
    # Doing b-score if version 2
    if (version == 2) {
      ntest <- ceiling(length(group_s) * borderp)
      scs <- numeric(ntest)
      for (j in 2:(ntest + 1)) {
        r1 <- sort(rand_rscores, decreasing = TRUE)[j-1]
        r2 <- sort(rand_rscores, decreasing = TRUE)[j]
        logp <- (N - length(group_s) + j) * (log(1 - r1) - log(1 - r2))
        scs[j-1] <- 1 - ((1 - r1) / (1 - r2))^(N - length(group_s) + j)
      }
      cs[counter] <- min(scs) * ntest
    }
  }
  
  if (return_mean) {
    return(mean(cs))
  } else {
    return(median(cs))
  }
  
}


