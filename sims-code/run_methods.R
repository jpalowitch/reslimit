source("recursive_mod.R")
sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
source("sims-code/netfuns.R")

set.seed(12345)
  
for (exper in exper_names) {

  # twocliq with increasing nS
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  for (p in seq_along(ps)) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
  
    for (i in 1:nsims) {
      
      # Making and saving data
      dat_fn <- paste0(i, ".dat")
      timer <- proc.time()[3]
      results <- recursive_mod(file.path(curr_dir, dat_fn))
      timer <- proc.time()[3] - timer
      save(results, timer, file = file.path(curr_dir, paste0("rmod_results_", i, ".RData")))
      
    }
    
  }
  
}
