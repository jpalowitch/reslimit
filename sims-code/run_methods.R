Args = commandArgs(trailingOnly=TRUE)

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

to_run <- as.numeric(Args)

source("recursive_mod.R")
source("cluster_resolution.R")
source("read_undirected_py.R")
source("sims-code/netfuns.R")
source("osRead.R")
do_rmod_louvain <- TRUE
do_rmod_simann <- FALSE
do_simann <- FALSE
sa_path <- "~/Documents/code/reslimit"
set.seed(12345)
  
for (exper in exper_names[to_run]) {

  # Loading/setting parameters
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  for (p in seq_along(ps)) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
  
    for (i in 1:nsims) {
      
      # Naming the files
      seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
      dat_fn <- paste0(i, ".dat")
      g_fn <- file.path(curr_dir, dat_fn)
      
      if (do_rmod_simann) {
        
        # Running recursive mod
        set.seed(as.numeric(readLines(seedfile)) + 1)
        timer <- proc.time()[3]
        res <- recursive_mod(g_fn, mod_type = "sim_ann", cscore_type = "r_cscore")
        timer <- proc.time()[3] - timer
        results <- res$results
        save(results, timer, file = file.path(curr_dir, paste0("rmodc_results_", i, ".RData")))
        
        # Running recursive mod w/bscore
        set.seed(as.numeric(readLines(seedfile)) + 2)
        timer <- proc.time()[3]
        res <- recursive_mod(g_fn, mod_type = "sim_ann", cscore_type = "r_bscore")
        timer <- proc.time()[3] - timer
        results <- res$results
        save(results, timer, file = file.path(curr_dir, paste0("rmodb_results_", i, ".RData")))
        
      }

      if (do_rmod_louvain) {
        
        # Running recursive mod with Louvain
        set.seed(as.numeric(readLines(seedfile)) + 1)
        timer <- proc.time()[3]
        res <- recursive_mod(g_fn, mod_type = "louvain", cscore_type = "r_cscore")
        timer <- proc.time()[3] - timer
        results <- res$results
        save(results, timer, file = file.path(curr_dir, paste0("rmodcLouv_results_", i, ".RData")))
        
        # Running recursive mod w/bscore with Louvain
        set.seed(as.numeric(readLines(seedfile)) + 1)
        timer <- proc.time()[3]
        res <- recursive_mod(g_fn, mod_type = "louvain", cscore_type = "r_bscore")
        timer <- proc.time()[3] - timer
        results <- res$results
        save(results, timer, file = file.path(curr_dir, paste0("rmodbLouv_results_", i, ".RData")))
        
      }
      
      # Running Louvain
      set.seed(as.numeric(readLines(seedfile)) + 3)
      timer <- proc.time()[3]
      G <- graph.edgelist(as.matrix(read.table(g_fn)), directed = FALSE)
      timer <- proc.time()[3] - timer
      comms <- cluster_louvain(G)$membership
      results <- lapply(1:max(comms), function (i) which(comms == i))
      save(results, timer, file = file.path(curr_dir, paste0("louvain_results_", i, ".RData")))
      
      
      if (do_simann) {
        
        # Running Simulated Annealing
        set.seed(as.numeric(readLines(seedfile)) + 4)
        timer <- proc.time()[3]
        savedir <- file.path(tempdir(), "sa_test")
        if (dir.exists(savedir)) unlink(savedir, recursive = TRUE)
        system(paste("python", file.path(sa_path, "clustering_programs_5_2/select.py"), 
                     "-n", file.path(curr_dir, dat_fn), "-p 8", 
                     paste("-f", savedir), "-c 1"))
        sa_results <- readLines(file.path(savedir, "results_1/tp"))
        sa_results <- lapply(sa_results, function (C) strsplit(C, "\t")[[1]])
        results <- lapply(sa_results, function (C) as.numeric(C))
        timer <- proc.time()[3] - timer
        save(results, timer, file = file.path(curr_dir, paste0("simann_results_", i, ".RData")))
        
      }
      
      # Running OSLOM
      set.seed(as.numeric(readLines(seedfile)) + 5)
      timer <- proc.time()[3]
      oldwd <- setwd(curr_dir)
      system2("./../../../OSLOM2/oslom_undir",
              c("-uw", "-f", dat_fn, "-fast", "-singlet"))
      setwd(oldwd)
      timer <- proc.time()[3] - timer
      results <- osRead(file.path(curr_dir, paste0(dat_fn, "_oslo_files"), "tp"))$communities
      save(results, timer, file = file.path(curr_dir, paste0("oslom_results_", i, ".RData")))
      
      # Running RB resolution tuning
      set.seed(as.numeric(readLines(seedfile)) + 6)
      timer <- proc.time()[3]
      res <- cluster_resolution(g_fn, res_start = 0, res_end = 10, interval = 0.1)
      timer <- proc.time()[3]
      
      

    }
    
  }
  
}
