source("recursive_mod.R")
sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
source("sims-code/netfuns.R")
source("osRead.R")

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
      
      # Naming the file
      if (truth_type == "LFR") {
        dat_fn <- paste0("network_", i, ".dat")
      } else {
        dat_fn <- paste0(i, ".dat")
      }

      # Running recursive mod
      timer <- proc.time()[3]
      results <- recursive_mod(file.path(curr_dir, dat_fn))
      timer <- proc.time()[3] - timer
      if (length(results[[1]]) > 0) {
        results <- lapply(1:max(results[[1]]), function (i) which(results[[1]] == i))
      } else {
        results <- NULL
      }
      save(results, timer, file = file.path(curr_dir, paste0("rmod_results_", i, ".RData")))
      
      # Running Louvain
      timer <- proc.time()[3]
      G <- graph.edgelist(as.matrix(read.table(file.path(curr_dir, dat_fn))), 
                          directed = FALSE)
      timer <- proc.time()[3] - timer
      comms <- cluster_louvain(G)$membership
      results <- lapply(1:max(comms), function (i) which(comms == i))
      save(results, timer, file = file.path(curr_dir, paste0("louvain_results_", i, ".RData")))
      
      # Running OSLOM
      timer <- proc.time()[3]
      oldwd <- setwd(curr_dir)
      system2("./../../../OSLOM2/oslom_undir",
              c("-uw", "-f", dat_fn, "-fast", "-singlet"))
      setwd(oldwd)
      timer <- proc.time()[3] - timer
      results <- osRead(file.path(curr_dir, paste0(dat_fn, "_oslo_files"), "tp"))$communities
      save(results, timer, file = file.path(curr_dir, paste0("oslom_results_", i, ".RData")))

    }
    
  }
  
}
