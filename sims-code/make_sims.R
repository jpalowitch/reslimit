sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
source("sims-code/netfuns.R")

# Remake sims or just write scripts?
remake_sims <- FALSE
oslom_exper_lines <- NULL

set.seed(12345)

for (exper in exper_names) {
  
  rootdir <- file.path(sim_res_dir, exper)
  if (!dir.exists(rootdir))
    dir.create(rootdir)
  
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  oslom_lines <- NULL
  
  for (p in 1:length(ps)) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
    oslom_lines <- c(oslom_lines, paste("cd", curr_dir))
    
    for (i in 1:nsims) {
      
      cat("p =", p, "i =", i, "\n")
      dat_fn <- paste0(i, ".dat")
      
      if (remake_sims) {
      
        # Setting seed
        seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
        if (!file.exists(seedfile)) {
          seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
          writeLines(seedpi, con = seedfile)
        }
        
        # Making and saving data
        if (make_type == "R_igraph") {
          set.seed(as.numeric(readLines(seedfile)))
          Gp <- eval(parse(text = make_code))
          write.table(get.edgelist(Gp), sep = "\t",
                      file = file.path(curr_dir, dat_fn),
                      quote = FALSE, row.names = FALSE, col.names = FALSE)
          save(Gp, file = file.path(curr_dir, paste0(i, ".RData")))
        }
        if (make_type == "system") {
          oldwd <- setwd(curr_dir)
          eval(parse(text = make_code))
          setwd(oldwd)
        }
      
        # Making truth
        if (truth_type == "Manual") {
          comms <- eval(parse(text = truth_code))
          truth_strings <- unlist(lapply(comms, paste, collapse = " "))
          writeLines(truth_strings,
                     con = file.path(curr_dir, paste0(i, "_truth.dat")))
        }
        if (truth_type == "LFR") {
          
          
        }
      
      }
      
      # Writing OSLOM lines
      oslom_add <- paste("./../../../OSLOM2/oslom_undir -uw -f", 
                         dat_fn, "-singlet", "-fast")
      oslom_lines <- c(oslom_lines, oslom_add)
      
    }
    
    oslom_lines <- c(oslom_lines, paste("cd", "../../../"))
    
  }

  writeLines(oslom_lines, con = file.path(rootdir, "oslom_lines.txt"))
  oslom_exper_lines <- c(oslom_exper_lines, paste("bash", file.path(rootdir, "oslom_lines.txt")))
  
}

writeLines(oslom_exper_lines, con = file.path(sim_res_dir, "oslom_exper_lines.txt"))
