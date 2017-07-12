recursive_mod <- function (fn, alpha = 0.05,
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
    cat("modularity maximization for level", level, "subgraph", sgn, "\n")
    G <- graph.edgelist(as.matrix(read.table(graph_fn)), directed = FALSE)
    n <- length(V(G))
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
    
    if (length(comms) > 1 && length(comms) < n) {
    
      # Compute c-b-scores
      signi_dir <- "signi-files"
      dir.create(signi_dir)
      setwd(signi_dir)
      system2(file.path(oldwd, compare_path, "compare"), 
              paste(paste("-f", file.path("..", graph_fn)), 
                    paste("-c", file.path("..", comms_fn)), 
                    "-t 0.01"))
      setwd("../")
      
    } else {
      
      writeLines(paste(1, n, 1, 1, 0, 0), con = paste0(sg_fn_base, ".dat.table"))
      
    }
    
    # Assess comms and save (if needed)
    cbscores <- read.table(paste0(sg_fn_base, ".dat.table"), sep = "", 
                           header = FALSE)
    which_write <- which(cbscores$V3 <= alpha)
    next_level_dir <- file.path("..", paste0("level", level + 1))
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
    
    # Assess background (if non-trivial) and save
    background <- unlist(comms[-which_write])
    if (length(background) > 0 && length(background) < n) {
      
      # Getting significance of background
      writeLines(paste(background, collapse = " "), con = "background.dat")
      setwd(signi_dir)
      system2(file.path(oldwd, compare_path, "compare"), 
              paste(paste("-f", file.path("..", graph_fn)), 
                    paste("-c", file.path("..", "background.dat")), 
                    "-t 0.01 -nobcore"))
      setwd("../")
      cbscores <- read.table(paste0(sg_fn_base, ".dat.table"), sep = "", 
                           header = FALSE)
      
      # If significant,
      if (cbscores$V3 <= 0.05) {
        
        # Add to comms and save
        comms <- c(comms, list(background))
        comms_text <- unlist(lapply(comms, paste, collapse = " "))
        writeLines(comms_text, con = comms_fn)
        which_write <- c(which_write, length(comms))
        
        # re-compute c-b-scores for future reading
        signi_dir <- paste0(sg_fn_base, "_signi")
        dir.create(signi_dir)
        setwd(signi_dir)
        system2(file.path(oldwd, compare_path, "compare"), 
                paste(paste("-f", file.path("..", graph_fn)), 
                      paste("-c", file.path("..", comms_fn)), 
                      "-t 0.01 -nobcore"))
        setwd("../")
        
        # write the background subgraph
        cur_next_nsgs <- length(list.files(next_level_dir))
        next_sg_base <- paste0("subgraph", cur_next_nsgs + 1)
        Gi <- delete.vertices(G, setdiff(1:n, comms[[length(comms)]]))
        write.table(get.edgelist(Gi), sep = "\t",
                    file = file.path(next_level_dir, 
                                     paste0(next_sg_base, ".dat")),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      }
        
    }
    
    # Creating membership vector
    membership <- numeric(n)
    for (i in seq_along(which_write)) {membership[comms[which_write][[i]]] <- i}
    return(membership)
    
  }
  
  # Level loop
  level <- 1
  level_memships <- list(NULL)
  repeat {
    
    # Prepping level
    level_dir <- file.path(mod_dir, paste0("level", level))
    if (level == 1) {
      dir.create(level_dir)
      file.copy(fn, to = level_dir)
      file.rename(list.files(level_dir, full.names = TRUE), file.path(level_dir, "subgraph1.dat"))
    }
    oldwd <- setwd(level_dir)
    
    # Running on subgraphs and 
    nsgs <- length(list.files())
    memships <- lapply(1:nsgs, level_run, level = level)
    
    # Check if next level was generated
    setwd(oldwd)
    if (!dir.exists(file.path(mod_dir, paste0("level", level + 1))))
      break
    
    # Saving memberships
    if (level == 1) {
      level_memships[[level]] <- unlist(memships)
    } else {
      membership <- numeric(length(level_memships[[1]]))
      for (j in seq_along(memships)) {
        cindx <- which(level_memships[[level - 1]] == j)
        membership[cindx] <- memships[[j]]
        nonzeros <- cindx[which(memships[[j]] != 0)]
        offset <- if (j > 1) {length(memships[[j - 1]])} else {0}
        membership[nonzeros] <- membership[nonzeros] + offset
      }
      level_memships[[level]] <- membership
      rm(membership)
    }
    
    level <- level + 1

  }
  
  return(level_memships)

}
