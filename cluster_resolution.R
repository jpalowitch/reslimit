cluster_resolution <- function (fn, method = 'RBConfiguration', res_start = 1, res_end = 1,
                                interval = (res_end - res_start) / 2, VI_thres = 0.005,
                                test_partitions = NULL) {
  
  # Loading libraries
  require(rPython)
  python.exec("import numpy")
  python.exec("import louvain")
  python.exec("import igraph as ig")
  
  # Getting graph
  python.exec(paste0("g = ig.Graph.Read_Ncol('", fn, "', directed = False)"))
  python.exec("n = g.vcount()")
  N <- python.get("n")
  
  if (res_start != res_end) {
    

    
    # Finding Ks and VIs
    res_pars <- seq(res_start, res_end, interval)
    VIs <- numeric(length(res_pars))
    Ks <- integer(length(res_pars))
    past_membership <- rep(1, N)
    
    # Setting up compare partitions, if any
    nparts <- length(test_partitions)
    if (nparts > 0) {
      part_mat <- matrix(0, nparts, length(res_pars))
    } else {
      part_mat <- NULL
    }
     
    # Res parameter loop   
    cat("going through res parameters...\n")
    for (i in seq_along(res_pars)) {
  
      python.exec(paste0("part = louvain.find_partition(g, method = '", method,
                         "', resolution_parameter = ", res_pars[i], ")"))
      python.exec("membership = part.membership")
      cur_membership <- python.get("membership") + 1
      
      if (nparts > 0) {
        for (j in 1:nparts) {
          testpart <- test_partitions[[j]]
          background <- which(testpart == 0)
          cur_membership_test <- cur_membership[-background]
          testpart <- testpart[-background]
          part_mat[j, i] <- compare(testpart, cur_membership_test, 
                                    method = "vi") / log(N)
        }
      }
      
      Ks[i] <- max(cur_membership)
      VIs[i] <- compare(past_membership, cur_membership, method = "vi") / log(N)
      past_membership <- cur_membership
      
    }
  
    # Finding range for Ks
    Krle <- rle(Ks)
    best_nc <- Krle$values[which.max(Krle$lengths)]
    K_locs <- which(Ks == best_nc)
    
    # Finding range for VIs
    VIs[VIs <= VI_thres] <- 0
    VI_subset <- VIs[K_locs]
    Vrle <- rle(VI_subset)
    zero_strings <- which(Vrle$values == 0)
    max_zero_string <- zero_strings[which.max(Vrle$lengths[zero_strings])]
    
    res_par_final <- res_pars[K_locs[max_zero_string]]
  
  } else {
    
    res_par_final <- res_start
    VIs <- NULL
    Ks <- NULL
    
  }
  
  # Finding final resolution parameter and returning
  python.exec(paste0("part = louvain.find_partition(g, method = '", method,
                     "', resolution_parameter = ", res_par_final, ")"))
  python.exec("membership = part.membership")
  membership <- python.get("membership") + 1
  
  # Returning
  return(list(membership = membership,
              res_par = res_par_final,
              part_mat = part_mat,
              res_pars = seq(res_start, res_end, interval),
              VIs = VIs, Ks = Ks))
  
}
