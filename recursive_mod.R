recursive_mod <- function (fn, alpha = 0.05, data_dir = ".",
                           compare_path = "cb-signi", 
                           mod_type = "louvain") {

  library(igraph)
  
  # Getting path to the file and making output directory
  fn_paths <- strsplit(fn, "/", fixed = TRUE)[[1]]
  if (length(fn_paths) > 1) {
    fn_dir <- paste(fn_paths[1:(length(fn_paths) - 1)], collapse = "/")
  } else {
    fn_dir <- "."
  }
  fn_base <- strsplit(tail(fn_paths, 1), ".", fixed = TRUE)[[1]][1]
  mod_dir <- file.path(fn_dir, paste0(fn_base, "_rmod_out"))
  if (dir.exists(mod_dir))
    unlink(mod_dir, recursive = TRUE)
  dir.create(mod_dir)
  
  # Function to analyze subgraph
  level_run <- function (sgn, level) {
    
    # Run modularity and save community file
    sg_fn_base <- paste0("subgraph", sgn)
    graph_fn <- file.path(paste0(sg_fn_base, ".dat"))
    G <- graph.edgelist(as.matrix(read.table(graph_fn)), directed = FALSE)
    if (mod_type == "louvain") {
      res <- cluster_louvain(G)
    } else {
      stop("mod_type unsupported\n")
    }
    comms <- lapply(1:max(res$membership), 
                    function (i) which(res$membership == i))
    comms_fn <- file.path(paste0(sg_fn_base, "_comms.dat"))
    comms_text <- unlist(lapply(comms, paste, collapse = " "))
    writeLines(comms_text, con = comms_fn)
    
    if (length(comms) > 1 && length(comms) < length(V(G))) {
    
      # Compute c-b-scores
      signi_dir <- paste0(sg_fn_base, "_signi")
      dir.create(signi_dir)
      setwd(signi_dir)
      system2(file.path(oldwd, compare_path, "compare"), 
              paste(paste("-f", file.path("../", graph_fn)), 
                    paste("-c", file.path("../", comms_fn)), 
                    "-t 0.01"))
      setwd("../")
      
    } else {
      
      writeLines(paste(1, length(V(G)), 1, 1, 0, 0), con = paste0(sg_fn_base, ".dat.table"))
      
    }
    
    # Assess comms and save (if needed)
    cbscores <- read.table(paste0(sg_fn_base, ".dat.table"), sep = "", 
                           header = FALSE)
    which_write <- which(cbscores$V4 <= alpha)
    next_level_dir <- file.path("../", paste0("level", level + 1))
    if (!dir.exists(next_level_dir) && length(which_write) > 0)
      dir.create(next_level_dir)
    cur_next_nsgs <- length(list.files(next_level_dir))
    for (i in seq_along(which_write)) {
      next_sg_base <- paste0("subgraph", cur_next_nsgs + i)
      Gi <- delete.vertices(G, setdiff(V(G), comms[[which_write[i]]]))
      write.table(get.edgelist(Gi), sep = "\t",
                  file = file.path(next_level_dir, paste0(next_sg_base, ".dat")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
  }
  
  # Level loop
  level <- 1
  
  repeat {
    
    # Prepping level
    level_dir <- file.path(mod_dir, paste0("level", level))
    if (level == 1) {
      dir.create(level_dir)
      file.copy(fn, to = level_dir)
      file.rename(list.files(level_dir, full.names = TRUE), file.path(level_dir, "subgraph1.dat"))
    }
    oldwd <- setwd(level_dir)
    
    # Running on subgraphs
    nsgs <- length(list.files())
    sapply(1:nsgs, level_run, level = level)
    
    # Check if next level was generated
    setwd(oldwd)
    level <- level + 1
    if (!dir.exists(file.path(mod_dir, paste0("level", level))))
      break
    
  }

}
