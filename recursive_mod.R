mycscore_path <- "~/Documents/code/reslimit"
source(file.path(mycscore_path, "my_cscore.R"))

recursive_mod <- function (fn, alpha = 0.05, cscore_type = "default",
                           compare_path = "cb-signi", sa_path = "~/Documents/code/reslimit",
                           mod_type = "louvain", borderp = 0.25) {

  library(igraph)
  
  # Getting path to the file and making output directory
  fn_paths <- strsplit(fn, "/", fixed = TRUE)[[1]]
  if (length(fn_paths) > 1) {
    fn_dir <- paste(fn_paths[1:(length(fn_paths) - 1)], collapse = "/")
  } else {
    fn_dir <- "."
  }
  fn_base <- strsplit(tail(fn_paths, 1), ".", fixed = TRUE)[[1]][1]
  mod_dir <- file.path(tempdir(), paste0(fn_base, "_rmod_out"))
  if (dir.exists(mod_dir))
    unlink(mod_dir, recursive = TRUE)
  dir.create(mod_dir)
  
################################################################################
  # Function to analyze subgraph
  level_run <- function (sgn, level) {
    
    # Loading network file and creating graph object
    sg_fn_base <- paste0("subgraph", sgn)
    graph_fn <- file.path(paste0(sg_fn_base, ".dat"))
    cat("modularity maximization for level", level, "subgraph", sgn, "\n")
    el <- as.matrix(read.table(graph_fn))
    G <- graph.edgelist(el[ , 1:2], directed = FALSE)
    n <- length(V(G))
    #cat("n =", n, "\n")
    
    # Add edges
    if (ncol(el) > 2) {
      E(G)$weight <- el[ , 3]
    } else {
      E(G)$weight <- rep(1, nrow(el))
    }
    
    # Running optimization
    if (mod_type == "louvain") {
      res <- cluster_louvain(G)
      comms <- lapply(1:max(res$membership), 
                    function (i) which(res$membership == i))
    } else if (mod_type == "spinglass") {
      res <- cluster_spinglass(G)
      comms <- lapply(1:max(res$membership), 
                    function (i) which(res$membership == i))
    } else if (mod_type == "sim_ann") {
      savedir <- "sa_test"
      if (dir.exists(savedir)) unlink(savedir, recursive = TRUE)
      system(paste("python", file.path(sa_path, "clustering_programs_5_2/select.py"), 
                   "-n", graph_fn, "-p 8", 
                   paste("-f", savedir), "-c 1"))
      sa_results <- readLines("sa_test/results_1/tp")
      sa_results <- lapply(sa_results, function (C) strsplit(C, "\t")[[1]])
      comms <- lapply(sa_results, function (C) as.numeric(C))
    } else {
      stop("Mod type unsupported")
    }
    #cat("n =", n, "\n")
    # Writing communities
    comms_fn <- file.path(paste0(sg_fn_base, "_comms.dat"))
    comms_text <- unlist(lapply(comms, paste, collapse = " "))
    writeLines(comms_text, con = comms_fn)
    
    # Getting c-scores
    if (length(comms) > 1 && length(comms) < n) {
    
      if (cscore_type %in% c("default", "bscore")) {
        
        # Compute original c/b-scores
        signi_dir <- "signi-files"
        if (!dir.exists(signi_dir)) dir.create(signi_dir)
        setwd(signi_dir)
        system2(file.path(oldwd, compare_path, "compare"), 
                paste(paste("-f", file.path("..", graph_fn)), 
                      paste("-c", file.path("..", comms_fn)), 
                      "-t 0.01"))
        setwd("../")
        cbscores <- read.table(paste0(sg_fn_base, ".dat.table"), sep = "", 
                               header = FALSE)
        
        # Need to input scores according to V1, so make a dummy vector
        cscores <- rep(1, length(comms))
        if (cscore_type == "default") {
          cscores[cbscores$V1] <- cbscores$V3
        } else {
          cscores[cbscores$V1] <- cbscores$V4
        }
      } else if(cscore_type %in% c("r_cscore", "r_bscore")) {
        
        # Compute my c/b scores
        if (cscore_type == "r_cscore") {
          cscores <- unlist(lapply(comms, my_cscore, G, 2))
        } else {
          cscores <- unlist(lapply(comms, my_cscore, G, 3, borderp))
        }
        
      } else {
        stop("cscore type unsupported\n")
      }
      
    } else {
      
      cscores <- 1
      writeLines(paste(1, n, 1, 1, 0, 0), con = paste0(sg_fn_base, ".dat.table"))
      
    }
    
    # Write cscores
    writeLines(as.character(cscores), con = "comm_cscores.txt")
    which_write <- which(cscores <= alpha)
    
    # Preparing next level (if needed)
    next_level_dir <- file.path("..", paste0("level", level + 1))
    if (!dir.exists(next_level_dir) && length(which_write) > 0)
      dir.create(next_level_dir)
    
    # Preparing subgraph writing and pointering
    cur_next_nsgs <- length(list.files(next_level_dir))
    comm_nums <- numeric(length(which_write))
    
    # The subgraph writing loop (note, will not execute if none significant)
    for (i in seq_along(which_write)) {
      
      # Naming the subgraph file
      next_sg_base <- paste0("subgraph", cur_next_nsgs + i)
      
      # Storing the community number to match the subgraph file
      comm_nums[i] <- cur_next_nsgs + i
      
      # Making the subgraph
      Gi <- delete.vertices(G, setdiff(V(G), comms[[which_write[i]]]))
      
      # Writing the subgraph
      write.table(get.edgelist(Gi), sep = "\t",
                  file = file.path(next_level_dir, paste0(next_sg_base, ".dat")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    # Making list of significant comms
    sig_comms <- comms[which_write]
    
    # Assess background (if non-trivial) and save
    background <- unlist(comms[setdiff(1:length(comms), which_write)])
    if (length(background) > 2 && length(background) < n) {
      
      
      if (cscore_type %in% c("default", "bscore")) {
        
        # Compute original c/b-scores
        writeLines(paste(background, collapse = " "), con = "background.dat")
        setwd(signi_dir)
        system2(file.path(oldwd, compare_path, "compare"), 
              paste(paste("-f", file.path("..", graph_fn)), 
                    paste("-c", file.path("..", "background.dat")), 
                    "-t 0.01"))
        setwd("../")
        bg_cbscores <- read.table(paste0(sg_fn_base, ".dat.table"), sep = "", 
                                  header = FALSE)
        
        if (cscore_type == "default") {
          bg_cscore <- bg_cbscores$V3
        } else {
          bg_cscore <- bg_cbscores$V4
        }
      } else if(cscore_type %in% c("r_cscore", "r_bscore")) {
        
        # Compute my c/b scores
        if (cscore_type == "r_cscore") {
          bg_cscore <- my_cscore(background, G, 2)
        } else {
          bg_cscore <- my_cscore(background, G, 3, borderp)
        }
        
      } else {
        stop("cscore type unsupported\n")
      }
      
      # If significant,
      if (bg_cscore <= 0.05) {
        
        # write the background subgraph
        cur_next_nsgs <- length(list.files(next_level_dir))
        next_sg_base <- paste0("subgraph", cur_next_nsgs + 1)
        comm_nums <- c(comm_nums, cur_next_nsgs + 1)
        Gi <- delete.vertices(G, unlist(comms[which_write]))
        write.table(get.edgelist(Gi), sep = "\t",
                    file = file.path(next_level_dir, 
                                     paste0(next_sg_base, ".dat")),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        # add to sig_comms
        sig_comms <- c(sig_comms, list(background))
        
      }
        
    }
    
    
    # Creating membership vector
    membership <- numeric(n)
    for (i in seq_along(sig_comms)) {membership[sig_comms[[i]]] <- comm_nums[i]}
    
    #cat("-- done with modularity maximization for level", level, "subgraph", sgn, "\n")
    return(membership)
    
  }
  
################################################################################
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
    
    # Running on subgraphs
    nsgs <- length(list.files())
    memships <- lapply(1:nsgs, level_run, level = level)
    
    # Check if next level was generated
    setwd(oldwd)
    if (!dir.exists(file.path(mod_dir, paste0("level", level + 1))))
      break
    
    # Saving memberships
    if (level == 1) {
      
      # Then just store the modularity maximization
      level_memships[[level]] <- unlist(memships)
      
    } else {
      
      # Setting up membership vector
      memship_vec <- numeric(length(level_memships[[1]]))
      
      for (j in seq_along(memships)) {
        
        # Finding out to which old community the new labels correspond
        cindx <- which(level_memships[[level - 1]] == j)
        
        # Putting the new labels in their position
        memship_vec[cindx] <- memships[[j]]
        
      }
      level_memships[[level]] <- memship_vec
      rm(memship_vec)
    }
    
    level <- level + 1

  }
  
  n <- length(memships[[1]])
    
  # Creating results list
  if (length(level_memships[[1]]) > 0) {
    
    # Collecting bottom-level result
  
    if (length(level_memships) > 1) {
      
      full_level_memships <- level_memships
      
      for (j in 2:length(level_memships)) {
      
        clust <- level_memships[[1]]
      
        for (lvl in 1:(j - 1)) {
        
          clust_start <- clust
          
          for (ci in sort(unique(clust_start))) {
            
            ci_locs <- which(clust_start == ci)
            ci_nums <- as.integer(level_memships[[lvl + 1]][ci_locs])
            next_locs <- which(clust_start > ci)
            if (sum(ci_nums != 0) > 0) {
              comm_nums <- sort(unique(ci_nums))
              comm_labels <- rank(comm_nums)
              comm_label_mch <- comm_labels[match(ci_nums, comm_nums)]
              clust[ci_locs] <- clust[ci_locs] + comm_label_mch - 1
              clust[next_locs] <- clust[next_locs] + max(comm_label_mch) - 1
            }
            
          }
        
        }
        
        full_level_memships[[j]] <- clust
        
      }
      
    } else {
      
      clust <- level_memships[[1]]
      full_level_memships <- list(clust)
      
    }
    
    results <- lapply(sort(unique(clust)), function (i) which(clust == i))
    background <- integer(0)
    if (sum(clust == 0) > 1) {
      background <- results[[1]]
      results <- results[-1]
    }
    
  } else {
    full_level_memships <- memships
    results <- NULL
    background <- 1:n
  }
  
  return(list(alllevels = full_level_memships,
              results = results,
              background = background))

}
