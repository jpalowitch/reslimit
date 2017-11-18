Args = commandArgs(trailingOnly=TRUE)

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

suppressWarnings(numeric_args <- as.numeric(Args))
endp <- NULL
startp <- 1
if (sum(is.na(numeric_args)) == 0) {
  to_run <- as.numeric(Args)
} else {
  startp_loc <- which(Args == "startp")
  startp <- numeric_args[startp_loc + 1]
  if (sum(Args == "endp") > 0) {
    endp_loc <- which(Args == "endp")
    endp <- numeric_args[endp_loc + 1]
  }
  to_run <- numeric_args[1:(startp_loc - 1)]
}

source("recursive_mod.R")
source("cluster_resolution.R")
source("read_undirected_py.R")
source("sims-code/netfuns.R")
source("osRead.R")
source("sbm_diagnostics.R")
do_rmod_louvain <- TRUE
do_rmod_simann <- FALSE
do_simann <- FALSE
sa_path <- "~/Documents/code/reslimit"
mod_bp_path <- "./modbp/mod"
set.seed(12345)
  
for (exper in exper_names[to_run]) {

  # Loading/setting parameters
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  # setting ps to run
  if (is.null(endp)) {
    endp <- length(ps)
  }
  
  for (p in startp:endp) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
  
    for (i in 1:nsims) {
      
      cat("doing exper", exper, "which is number", match(exper, exper_names), "p =", p, "i = ", i, "\n")
      
      # Naming the files
      seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
      dat_fn <- paste0(i, ".dat")
      g_fn <- file.path(curr_dir, dat_fn)
      
      #edgelist <- read.table(g_fn)
      #G <- graph.edgelist(as.matrix(edgelist), directed = FALSE)
      #adjMat <- get.adjacency(G)
      #plot(adjMat[nodeorder, nodeorder])
      
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
      results <- lapply(1:max(comms), function (j) which(comms == j))
      save(results, timer, file = file.path(curr_dir, paste0("louvain_results_", i, ".RData")))
      
      # Running Infomap
      set.seed(as.numeric(readLines(seedfile)) + 3)
      timer <- proc.time()[3]
      G <- graph.edgelist(as.matrix(read.table(g_fn)), directed = FALSE)
      timer <- proc.time()[3] - timer
      comms <- cluster_infomap(G)$membership
      results <- lapply(1:max(comms), function (j) which(comms == j))
      save(results, timer, file = file.path(curr_dir, paste0("infomap_results_", i, ".RData")))
      
      
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
      cr_results <- res$membership
      results <- lapply(1:max(cr_results), function (j) which(cr_results == j))
      res_par <- res$res_par
      save(results, res_par, timer, file = file.path(curr_dir, paste0("RBres_results_", i, ".RData")))
      
      # Running CPM resolution tuning
      set.seed(as.numeric(readLines(seedfile)) + 7)
      timer <- proc.time()[3]
      res <- cluster_resolution(g_fn, res_start = 0, res_end = 1, interval = 0.01, method = "CPM")
      timer <- proc.time()[3]
      cr_results <- res$membership
      results <- lapply(1:max(cr_results), function (j) which(cr_results == j))
      res_par <- res$res_par
      save(results, res_par, timer, file = file.path(curr_dir, paste0("CPMres_results_", i, ".RData")))
      
      # Running MODRB
      modrb_seed <- as.numeric(readLines(seedfile)) + 8
      timer <- proc.time()[3]
      system(paste(mod_bp_path, "hiera", "-l", file.path(curr_dir, paste0(i, ".gml")),
                   "--confi", file.path(curr_dir, "modbp_result.dat")))
      membership <- as.numeric(unlist(readLines(file.path(curr_dir, "modbp_result.dat"))))
      results <- lapply(1:max(membership), function (j) which(membership == j))
      timer <- proc.time()[3] - timer
      save(results, timer, file = file.path(curr_dir, paste0("modbpy_results_", i, ".RData")))
      
      # running graphtool
      timer <- proc.time()[3]
      system(paste("/usr/bin/python",
                   "fit_sbm.py",
                   file.path(curr_dir, paste0(i, ".graphml")),
                   file.path(curr_dir, paste0(i, "_gtMemship.dat"))))
      timer <- proc.time()[3] - timer
      membership <- as.integer(readLines(file.path(curr_dir, paste0(i, "_gtMemship.dat")))) + 1
      results <- lapply(1:max(membership), function (j) which(membership == j))
      save(results, timer, file = file.path(curr_dir, paste0("sbm_results_", i, ".RData")))
      
      # For SBMs, assessing graph statistics
      if (grepl('SBM', exper)) {
        truth_fn <- file.path(curr_dir, paste0(i, "_truth.dat"))
        results <- SBM_diagnostics(G, truth_fn)
        save(results, file = file.path(curr_dir, paste0("sbm_diagnostics_", i, ".RData")))
      }
      

    }
    
  }
  
}
