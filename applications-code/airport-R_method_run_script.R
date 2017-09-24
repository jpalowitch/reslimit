saveDir <- "applications-results/airports"
source("saveit.R")

origwd <- setwd("../CCME")
source("CCME.R")
setwd(origwd)

library(igraph)

monthNames <- c("jan", "feb", "mar", "apr", 
                "may", "jun", "jul", "aug", 
                "sep", "oct", "nov", "dec")

fnList <- c(monthNames, "year")
load("applications-results/airports/data/adjList_cleaned.RData")
startYear <- as.numeric(names(adjList)[1])
endYear <- as.numeric(names(adjList)[length(adjList)])
yL <- endYear - startYear + 1

# Initializing run lines for SLPAw
slpa_run_lines <- character(0)
oslom_run_lines <- character(0)
fullNames <- character(0)
datNames <- character(0)

# Tracker for debugging
onMethod <- NULL

# Making the OSLOM folder
oslomdir <- "applications-results/airports/OSLOM2"
if (!dir.exists(oslomdir))
  dir.create(oslomdir, recursive = TRUE)

set.seed(12345)

for (y in 1:yL) {
    
  yearNum <- startYear + y - 1
  curr_dir <- file.path(saveDir, yearNum)
  fnList_y <- c(monthNames[1:length(adjList[[y]]$months)], "year")
  
  for (fn in fnList_y) {
      
    load(file.path(curr_dir, paste0(fn, ".RData")))
  
    # Formatting data for CCME run
    colnames(data) <- NULL
    rownames(data) <- NULL
    G <- graph.adjacency(data, mode = "undirected", weighted = "weight")
    edge_list <- get.edgelist(G)
    edge_list <- cbind(edge_list, E(G)$weight)
    edge_list <- as.data.frame(edge_list)
    names(edge_list) <- c("node1", "node2", "weight")
    
    dat_fn <- paste0(yearNum, "_", fn, ".dat")
    
    write.table(edge_list,
                file = file.path(oslomdir, dat_fn),
                row.names = FALSE,
                col.names = FALSE)
    
    # Setting seed
    seedfn <- file.path(curr_dir, paste0("seed_ccme_", fn, ".txt"))
    if (file.exists(seedfn)) {
      seed_draw <- as.integer(readLines(seedfn))
    } else {
      seed_draw <- sample(1e6, 1)
      writeLines(as.character(seed_draw), con = seedfn)
    }
    set.seed(seed_draw)

    # Running  
    Timer <- proc.time()[3]
    results <- CCME(edge_list)
    Timer <- proc.time()[3] - Timer
    results$time <- Timer
    save(results, file = file.path(curr_dir, 
                                   paste0("output_ccme_",fn,".RData")))
    
    # Running  
    Timer <- proc.time()[3]
    results <- CCME(edge_list, fastInitial = TRUE)
    Timer <- proc.time()[3] - Timer
    results$time <- Timer
    save(results, file = file.path(curr_dir, 
                                   paste0("output_ccme_fast_",fn,".RData")))
    
    # Running  
    Timer <- proc.time()[3]
    results <- CCME(edge_list, fastInitial = FALSE, zoverlap = TRUE)
    Timer <- proc.time()[3] - Timer
    results$time <- Timer
    save(results, file = file.path(curr_dir, 
                                   paste0("output_ccme_z_",fn,".RData")))
  
    
    # Making an ipairs file and saving
    write.table(edge_list, 
                file = file.path(curr_dir, paste0("network_", fn, ".ipairs")), 
                row.names = FALSE, 
                col.names = FALSE)
    
    # Adding line to SLPA script
    seed_draw <- sample(1e6, 1)
    slpa_addline <- paste0("java -jar ",
                           file.path( "methodFiles/GANXiS_v3.0.2/GANXiSw.jar"),
                           " -i ",
                           file.path(curr_dir, 
                                     paste0("network_", fn, ".ipairs")),
                           " -Sym 1 -r 0.1 -d ",
                           file.path(curr_dir),
                           " -seed ",
                           seed_draw)
    slpa_run_lines <- c(slpa_run_lines,
                        slpa_addline)
    
    # Adding OSLOM line to script
    inside_addline <- paste0("./oslom_undir -f ",
                             dat_fn,
                             " -w -singlet")
    oslom_run_lines <- c(oslom_run_lines,
                         inside_addline)
    
    fullNames <- c(fullNames, file.path(curr_dir, fn))
    datNames <- c(datNames, dat_fn)
    
      
  }

}

runcodedir <- "applications-results/airports/run-code"

if (!dir.exists(runcodedir))
  dir.create(runcodedir, recursive = TRUE)

# Initializing scripts to run batches, and logfiles
oslomfn0 <- file.path(runcodedir, "run-oslom")
slpafn0 <- file.path(runcodedir, "run-slpa")
file.create(paste0(c(oslomfn0, slpafn0), ".txt"))
file.create(paste0(c(oslomfn0, slpafn0), "_log.txt"))

# Saving slpa lines
fileConn <- file(paste0(slpafn0, ".txt"), "a")
writeLines(slpa_run_lines, con = fileConn)
close(fileConn)

# Saving oslom lines, plus the lines to copy and compile the program
fileConn <- file(paste0(oslomfn0, ".txt"), "a")

    # Coping OSLOM files
    writeLines(paste0("cp -R methodFiles/OSLOM2/* ", oslomdir),
               fileConn)
    
    # Set the directory
    writeLines(paste0("cd ", oslomdir),
               fileConn)
    
    # Compiling OSLOM
    writeLines("chmod 744 compile_all.sh && ./compile_all.sh",
               fileConn)
    
    # Writing run lines
    writeLines(oslom_run_lines, con = fileConn)
      
    # Re-set directory
    writeLines(paste0("cd ../../../"),
               fileConn)
    
close(fileConn)

# Saving oslom filenames
save(datNames, fullNames, file = file.path(oslomdir, "aports.RData"))
