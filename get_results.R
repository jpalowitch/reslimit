source("recursive_mod.R")

set.seed(12345)

# twocliq with increasing nS
rootdir <- "sims/twocliq_nS"
ps <- round(50 * 3^(seq(0, 4, 0.5)))
nsims <- 3

timers <- numeric(length(ps))

for (p in 1:6) {
  
  # Setting directory
  curr_dir <- file.path(rootdir, p)
  if (!dir.exists(curr_dir))
    dir.create(curr_dir)

  for (i in 1:nsims) {
    
    # Setting seed
    seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
    if (!file.exists(seedfile)) {
      seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
      writeLines(seedpi, con = seedfile)
    }
    
    # Loading results
    load(file.path(curr_dir, "rmod_results.RData"))
    timers[p] <- timers[p] + timer / nsims
    
  }
  
}
