Args = commandArgs(trailingOnly=TRUE)

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

to_run <- as.numeric(Args)

source("sims-code/netfuns.R")
source("full_dcsbm.R")
source("igraph_recast.R")
source("generate_sbm.R")
# Remake sims or just write scripts?
remake_sims <- TRUE
set.seed(12345)

for (exper in exper_names[to_run]) {
  
  rootdir <- file.path(sim_res_dir, exper)
  if (dir.exists(rootdir)) unlink(rootdir, recursive = TRUE)
  dir.create(rootdir)
  
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  for (p in 1:length(ps)) {
    
    # Setting directory
    curr_dir <- file.path(rootdir, p)
    if (!dir.exists(curr_dir))
      dir.create(curr_dir)
    
    for (i in 1:nsims) {
      
      cat("exper =", exper, "\n")
      cat("-- p =", p, "i =", i, "\n")
      dat_fn <- paste0(i, ".dat")
      g_fn <- file.path(curr_dir, dat_fn)
      
      if (remake_sims) {
      
        # Getting seed
        seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
        if (!file.exists(seedfile)) {
          seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
          writeLines(seedpi, con = seedfile)
        }
        
        # Making and saving data
        if (make_type == "R_igraph") {
          set.seed(as.numeric(readLines(seedfile)))
          dummy <- sapply(make_code, function (c) eval(parse(text = c)))
          write.table(get.edgelist(Gp), sep = "\t", file = g_fn,
                      quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
        if (make_type == "System") {
          set.seed(as.numeric(readLines(seedfile)))
          oldwd <- setwd(curr_dir)
          dummy <- sapply(make_code, function (c) eval(parse(text = c)))
          setwd(oldwd)
        }
        
        # Making truth in list format
        if (truth_type == "Manual") {
          dummy <- sapply(truth_code, function (c) eval(parse(text = c)))
        }
        if (truth_type == "LFR") {
          comms <- read.table(file.path(curr_dir, "community.dat"),
                              header = FALSE, sep = "\t")
          comms <- lapply(1:max(comms$V2), function (j) {
            return(comms$V1[comms$V2 == j])
          })
        }
        
        # Making membership vector
        membership <- unlist(comms)
        for (k in 1:length(comms)) membership[comms[[k]]] <- k
        
        # Recasting in python-igraph format
        cast_obj <- igraph_recast(as.matrix(read.table(g_fn)), membership)
        write.table(cast_obj$edgelist, sep = "\t", file = g_fn,
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        # Porting to python-igraph indexing
        comms <- lapply(comms, function (C) cast_obj$lookup[C])
        
        # Writing in GML format
        G_igraph <- graph.edgelist(cast_obj$edgelist, directed = FALSE)
        write_graph(G_igraph, file = file.path(curr_dir, paste0(i, ".gml")),
                    format = "gml")
        
        # Writing in graphml format
        write_graph(G_igraph, file = file.path(curr_dir, paste0(i, ".graphml")),
                    format = "graphml")
        
        # Saving truth list elements in lines of a dat file
        truth_strings <- unlist(lapply(comms, paste, collapse = " "))
        writeLines(truth_strings,
                   con = file.path(curr_dir, paste0(i, "_truth.dat")))
      
      }
      
    }
    
  }
  
}
