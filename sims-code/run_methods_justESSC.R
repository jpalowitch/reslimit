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
library(ESSC)
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
      
      edgelist <- read.table(g_fn)
      G <- graph.edgelist(as.matrix(edgelist), directed = FALSE)
      adjMat <- get.adjacency(G)
      #plot(adjMat[nodeorder, nodeorder])
      
      # Running recursive mod w/bscore
      set.seed(as.numeric(readLines(seedfile)) + 10)
      timer <- proc.time()[3]
      res <- essc(adjMat, alpha = 0.05, Null = "Binomial")
      timer <- proc.time()[3] - timer
      results <- res$Communities
      if (is.null(results[[1]])) results <- NULL
      save(results, timer, file = file.path(curr_dir, paste0("essc_results_", i, ".RData")))
      
      truth <- readLines(file.path(curr_dir, paste0(i, "_truth.dat")))
      comms <- lapply(truth, function (L) as.numeric(strsplit(L, " ", fixed = TRUE)[[1]]))

    }
    
  }
  
}
