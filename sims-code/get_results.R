Args = commandArgs(trailingOnly=TRUE)

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

to_run <- as.numeric(Args)

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
  rm(dummy_dir, fns, methFiles)
  
  # Making results matrices
  timer_mat <- onmi_mat <- matrix(0, length(methNames), length(ps))
  rownames(timer_mat) <- rownames(onmi_mat) <- methNames
  
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
        
        # Loading results and computing onmi
        load(res_rfn)
        if (length(results) > 0) {
          comm_strings <- unlist(lapply(results, paste, collapse = " "))
          writeLines(comm_strings, res_dfn)
          system2(mutual3, c(res_dfn, truth_fn, paste0("> ", mut_dfn)))
          onmi <- as.numeric(strsplit(readLines(mut_dfn), "\t")[[1]][2])
        } else {
          writeLines("", res_dfn)
          onmi <- 0
        }
        
        # Storing results
        timer_mat[meth, p] <- timer_mat[meth, p] + timer / nsims
        onmi_mat [meth, p] <- onmi_mat [meth, p] + onmi  / nsims
        
      }
      
    }
    
  }
  
  write.table(timer_mat, quote = FALSE, col.names = FALSE,
              file = file.path(rootdir, "timer_mat.txt"))
  write.table(onmi_mat, quote = FALSE, col.names = FALSE,
  file = file.path(rootdir, "onmi_mat.txt"))
  
}
