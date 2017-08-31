sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
source("sims-code/netfuns.R")

# Remake sims or just write scripts?
remake_sims <- TRUE
set.seed(12345)

for (exper in exper_names[-1]) {
  
  rootdir <- file.path(sim_res_dir, exper)
  if (!dir.exists(rootdir))
    dir.create(rootdir)
  
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  for (p in 1:length(ps)) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
    
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
        if (make_type == "System") {
          oldwd <- setwd(curr_dir)
          sapply(make_code, function (c) eval(parse(text = c)))
          setwd(oldwd)
        }
      
        # Making truth in list format
        if (truth_type == "Manual") {
          comms <- eval(parse(text = truth_code))
        }
        if (truth_type == "LFR") {
          comms <- read.table(file.path(curr_dir, "community.dat"),
                              header = FALSE, sep = "\t")
          comms <- lapply(1:max(comms$V2), function (j) {
            return(comms$V1[comms$V2 == j])
          })
        }
        
        # Saving truth list elements in lines of a dat file
        truth_strings <- unlist(lapply(comms, paste, collapse = " "))
        writeLines(truth_strings,
                   con = file.path(curr_dir, paste0(i, "_truth.dat")))
      
      }
      
    }
    
  }
  
}
