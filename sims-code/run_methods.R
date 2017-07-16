source("recursive_mod.R")

set.seed(12345)

# twocliq with increasing nS
rootdir <- "sims/twocliq_nS"
ps <- round(50 * 3^(seq(0, 4, 0.5)))
nsims <- 3


for (p in 1:(length(ps) - 1)) {
  
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
    
    # Making and saving data
    set.seed(as.numeric(readLines(seedfile)))
    dat_fn <- paste0(i, ".dat")
    timer <- proc.time()[3]
    results <- recursive_mod(file.path(curr_dir, dat_fn))
    timer <- proc.time()[3] - timer
    save(results, timer, file = file.path(curr_dir, paste0("rmod_results_", i, ".RData")))
    
  }
  
}
