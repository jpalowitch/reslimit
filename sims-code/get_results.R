Args = commandArgs(trailingOnly=TRUE)

# Compute times in log values?
log_time <- TRUE

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

to_run <- as.numeric(Args)

# Methods to ignore
ignore_meths <- c("rmodcLouv", "CPMres")

mutual3 <- "./mutual3/mutual"

for (exper in exper_names[to_run]) {
  
  # Loading/setting parameters
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  # Finding methNames
  dummy_dir <- file.path(rootdir, 1)
  fns <- list.files(dummy_dir)
  methFiles <- fns[grepl("results", fns)]
  methNames <- unique(sapply(methFiles, function (c) strsplit(c, "_")[[1]][1]))
  methNames <- setdiff(methNames, ignore_meths)
  rm(dummy_dir, fns, methFiles)
  
  # Making results matrices
  timer_mat <- onmi_mat <- CIB_mat <- matrix(0, length(methNames), length(ps))
  rownames(timer_mat) <- rownames(onmi_mat) <- rownames(CIB_mat) <- methNames
  
  for (p in seq_along(ps)) {
    
    cat("experiment", exper, "p =", p, "\n")
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)

    for (meth in methNames) {
  
      for (i in 1:nsims) {
        
        truth_fn <- file.path(curr_dir, paste0(i, "_truth.dat"))
        res_rfn  <- file.path(curr_dir, paste0(meth, "_results_", i, ".RData"))
        res_dfn  <- file.path(curr_dir, paste0(meth, "_results_", i, ".dat"))
        mut_dfn  <- file.path(curr_dir, paste0(meth, "_onmi_", i, ".dat"))
        
        # Loading truth
        truth <- readLines(truth_fn)
        truth <- lapply(truth, function (C) as.numeric(strsplit(C, " ")[[1]]))
        N <- length(unlist(truth))
        
        
        
        # Loading results and adjusting time if needed
        load(res_rfn)
        if (log_time) timer <- log10(timer + 1)
        
        # Computing onmi
        if (length(results) > 0) {
          
          # Reducing truth
          all_nodes <- unlist(results)
          truth <- lapply(truth, function (C) intersect(C, all_nodes))
          truth <- truth[!(lapply(truth, length) == 0)]
          truth <- unlist(lapply(truth, function (C) as.character(paste(C, collapse = " "))))
          meth_truth_fn <- file.path(curr_dir, paste0(meth, "_truth_", i, ".dat"))
          
          # Writing
          writeLines(truth, con = meth_truth_fn)
          comm_strings <- unlist(lapply(results, paste, collapse = " "))
          writeLines(comm_strings, res_dfn)
          system2(mutual3, c(res_dfn, meth_truth_fn, paste0("> ", mut_dfn)))
          onmi <- as.numeric(strsplit(readLines(mut_dfn), "\t")[[1]][2])
          CIB <- (N - length(unique(unlist(results)))) / N
          
        } else {
          writeLines("", res_dfn)
          onmi <- 0
          CIB <- 1
        }
        
        # Storing results
        timer_mat[meth, p] <- timer_mat[meth, p] + timer / nsims
        onmi_mat [meth, p] <- onmi_mat [meth, p] + onmi  / nsims
        CIB_mat  [meth, p] <- CIB
        
      }
      
    }
    
  }
  
  write.table(timer_mat, quote = FALSE, col.names = FALSE,
              file = file.path(rootdir, "timer_mat.txt"))
  write.table(onmi_mat, quote = FALSE, col.names = FALSE,
              file = file.path(rootdir, "onmi_mat.txt"))
  write.table(CIB_mat, quote = FALSE, col.names = FALSE,
              file = file.path(rootdir, "CIB_mat.txt"))
  
}
