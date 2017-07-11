source("osRead.R")

rootdir <- file.path("sims", "twocliq_ns")
load(file.path(rootdir, "pars.RData"))

for (p in 1:length(ps)) {
  
  # Setting directory
  curr_dir <- file.path(rootdir, p)
  
  for (i in 1:nsims) {
    
    dat_fn <- paste0(i, ".dat")
    tp_fn <- file.path(curr_dir, paste0(dat_fn, "_oslo_files"), "tp")
    results <- osRead(tp_fn)
    res_fn <- file.path(curr_dir, paste0("oslom_results_", i, ".RData"))
    save(results, file = res_fn)
    
  }
  
}