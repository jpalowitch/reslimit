library(igraph)
library(gdata)
power_law <- function (n, tau, max_set, min_set = NULL, mean_set = NULL) {
  
  if (is.null(mean_set) && is.null(min_set)) {
    cat("err: one of min_set, mean_set must be specified\n")
    return(NULL)
  }
  
  if (is.null(mean_set)) {
    prob0 <- (min_set:max_set)^(-tau)
    prob <- prob0 / sum(prob0)
  }
  
  if (is.null(min_set)) {
    min_set <- max_set
    prob0 <- (min_set:max_set)^(-tau)
    prob <- prob0 / sum(prob0)
    expv <- sum(prob * min_set:max_set)
    while (expv > mean_set) {
      min_set <- min_set - 1
      prob0 <- (min_set:max_set)^(-tau)
      prob <- prob0 / sum(prob0)
      expv <- sum(prob * min_set:max_set)
    }
    if (min_set < 1) {
      cat("err: mean is too small\n")
      return(NULL)
    }
  }
  
  if (min_set == max_set) {
    p_law_sample <- rep(max_set, n)
  } else {
    p_law_sample <- sample(min_set:max_set, n, replace = TRUE, prob = prob)
  }
  return(p_law_sample)
}

# Function to make a paramList
make_param_list <- function (N = 5000,
                             max_c = N * 3 / 10,
                             min_c = max_c * 2 / 3,
                             k = round(sqrt(N)),
                             min_k = NULL,
                             max_k = 3 * k,
                             tau_1 = 2,
                             tau_2 = 1,
                             s2n = 3,
                             mu = 0.3) {
  return(list("N" = N,
              "max_c" = max_c,
              "min_c" = min_c,
              "k" = k,
              "min_k" = min_k,
              "max_k" = max_k,
              "tau_1" = tau_1,
              "tau_2" = tau_2,
              "s2n" = s2n,
              "mu" = mu))
}

DCSBM <- function (param_list, type = "fast",
                   degrees = NULL, P = NULL,
                   membership = NULL,
                   muversion = TRUE, adjust = TRUE,
                   graph_tool = FALSE) {
  
  # Making a notes character vector
  notes <- character(0)
  
  # Extracting from paramList  
  N <- ifelse(is.null(membership), param_list$N, length(membership))
  max_c <- param_list$max_c
  min_c <- param_list$min_c
  min_k <- param_list$min_k
  k <- ifelse(is.null(membership), param_list$k, nrow(P))
  max_k <- param_list$max_k
  tau_1 <- param_list$tau_1
  tau_2 <- param_list$tau_2
  s2n <- param_list$s2n
  mu <- param_list$mu
  
  #cat("got here 1\n")
  
  # Making membership
  if (is.null(membership)) {
    
    # Adding communities until amount of comms is >= o_m + 1
    # and until the sizes sum to greater than n_m
    n_m <- N
    comm_sizes <- integer(0)
    n_c <- 0
    N_comm <- sum(comm_sizes)
    while (n_c < 1 || N_comm < n_m) {
      
      if (tau_2 < Inf) {
        new_draw <- power_law(1,
                              max_set = max_c,
                              min_set = min_c,
                              tau = tau_2)
      } else {
        if (min_c < max_c) {
          new_draw <- sample(min_c:max_c, 1)
        } else {
          new_draw <- max_c
        }
      }
      comm_sizes <- c(new_draw, comm_sizes)
      n_c <- length(comm_sizes)
      N_comm <- sum(comm_sizes)
    }
    
    # Taking nodes away proportionally so that
    # n_m is approx. N_comm
    comm_props <- comm_sizes / N_comm
    diff <- N_comm - n_m
    comm_sizes <- comm_sizes - round(comm_props * diff)
    N_comm <- sum(comm_sizes)
    
    # If N_comm still != n_m, correct by largest community
    diff <- N_comm - n_m
    comm_sizes[which.max(comm_sizes)] <- max(comm_sizes) - diff
    
    # Give memberships to nodes
    n_memberships <- rep(1, N)
    expanded_node_list <- unlist(lapply(1:N, function (node) {
      rep(node, n_memberships[node])}))
    comm_indx_vec <- unlist(lapply(1:n_c, function (cindx) {
      rep(cindx, comm_sizes[cindx])}))
    comm_indx_vec <- sample(comm_indx_vec)
    get_membership <- function (node) {
      comm_indx_vec[expanded_node_list == node]
    }
    membership <- unlist(lapply(1:N, get_membership))
    
  }
  #cat(head(membership), "\n")
  #cat("got here 2\n")
  
  K <- max(membership)
  if (is.null(P)) {
    # Making P matrix
    P <- matrix(1, K, K)
    if (!muversion) {
      diag(P) <- s2n
    } else {
      perComm <- mu / (K - 1)
      diag(P) <- 1 - mu
      P[row(P) != col(P)] <- perComm
    }
    P <- K * P / sum(P)
    if (type == "slow") P <- P * K
  }
  
  #cat("got here 3\n")
  
  # Setting degrees
  if (is.null(degrees)) {
    if (!is.null(min_k)) {
      degrees <- power_law(N, tau_1, max_set = max_k, min_set = min_k)
    } else {
      degrees <- power_law(N, tau_1, max_set = max_k, mean_set = k)
    }
  } 
  dT <- sum(degrees)
  
  if (!graph_tool) {
  
    cat("!graph_tool\n")
    
    # Fixing degrees step 1
    if (max(degrees) > N) {
      message("adjusting degrees, since max(degrees) > N\n")
      degrees[degrees > N] <- N
    }
    #cat("got here 4\n")
     
    while (adjust) { 
      
      # Fixing degrees step 2
      {
        # Getting unique community labels
        clabs <- sort(unique(membership))
        
        # Looping degree correction until fixed
        max_Pnum <- 2
        while (max_Pnum > 1) {
          
          # Finding maximum degree in each community
          maxd_C <- unlist(lapply(clabs, function (L) {
            max(degrees[membership == L])
          }))
          
          # Finding 2nd max degree in each community
          maxd_C2 <- unlist(lapply(clabs, function (L) {
            tail(sort(degrees[membership == L]), 2)[1]
          }))
          
          # Making max probability mat
          max_P <- tcrossprod(maxd_C) / dT * P
          diag(max_P) <- diag(max_P) * maxd_C2 / maxd_C
          
          # Assessing the max
          max_Pnum <- max(max_P)
          
          #cat(max_Pnum, "\n")
          
          #cat("got here 4.5\n")
          
          # Fixing if needed
          if (max_Pnum > 1) {
            culprit <- which(max_P == max_Pnum, arr.ind = TRUE)
            if (length(culprit) == 2) {
              C1 <- culprit[1]; C2 <- culprit[2]
            } else {
              C1 <- culprit[1, 1]; C2 <- culprit[1, 2]
            }
            D1 <- degrees[membership == C1]; D2 <- degrees[membership == C2]
            D1max <- max(D1); D2max <- max(D2)
            if (D1max >= D2max) {
              D1[D1 == D1max] <- D1max - 1
              degrees[membership == C1] <- D1
            } else {
              D2[D2 == D2max] <- D2max - 1
              degrees[membership == C2] <- D2
            }
          }
            
        }
        
      }
      
      # Scaling P so that it matches dT
      target_sum <- sum(degrees)
      dC <- unlist(lapply(1:K, function (k) sum(degrees[membership == k])))
      mean_mat <- tcrossprod(dC) / dT * P
      actual_sum <- sum(mean_mat)
      P <- P * target_sum / actual_sum
      #cat("iteration--\n")
      #cat("--", round(target_sum, 2), "vs", round(actual_sum, 2), "\n")
      if (actual_sum - target_sum <= 1e-5) break
    }
    #cat("got here 5\n")
    # Making block model
    cat("constructing block model...\n")
    
    if (type == "fast") {
    
      # Find relative inter-community weights
      comm_list <- lapply(clabs, function (L) which(membership == L))
      dC <- unlist(lapply(comm_list, function (C) sum(degrees[C])))
      Dmat <- tcrossprod(dC) * P / dT
      Dmat <- round(Dmat / sum(Dmat) * dT)
      
      cat("sum(Dmat) =", sum(Dmat), "\n")
      
      # Make number of edges
      nEdges <- rbinom(1, N^2, sum(Dmat / N^2))
      
      cat("nEdges =", nEdges, "\n")
      
      # Assigning edges to inter-community relationships
      Pmat <- Dmat / sum(Dmat)
      upperTriangle(Pmat) <- upperTriangle(Pmat) * 2
      lowerTriangle(Pmat) <- 0
      edgeCounts <- rmultinom(1, nEdges, as.vector(Pmat))
      edgeCountMat <- matrix(edgeCounts, K, K)
      upperTriangle(edgeCountMat) <- upperTriangle(edgeCountMat) / 2
      lowerTriangle(edgeCountMat) <- lowerTriangle(t(edgeCountMat))
      edgeCountMat <- ceiling(edgeCountMat)
      
      # Assigning inter-community edge counts to nodes
      nodeDegMat <- matrix(0, N, K)
      for (i1 in 1:K) {
        nodes1 <- comm_list[[i1]]
        degs1 <- degrees[nodes1]
        for (i2 in i1:K) {
          nodeDegMat[nodes1, i2] <- rmultinom(1, edgeCountMat[i1, i2], degs1)
        }
      }
      
      # Run within/between-community configuration models
      edgelist_mats <- rep(list(NULL), K^2)
      for (i1 in 1:K) {
        
        nodes1 <- comm_list[[i1]]
        
        for (i2 in i1:K) {
          
          elm_indx <- K * (i1 - 1) + i2
          
          # Draw nodes node-wise for between-comm
          if (i2 > i1) {
            nodes2 <- comm_list[[i2]]
            node2_input <- unlist(lapply(nodeDegMat[nodes1, i2], function (d) {
              sample(nodes2, d, replace = FALSE, prob = degrees[nodes2])
            }))
            node1_input <- rep.int(nodes1, times = nodeDegMat[nodes1, i2])
            edgelist_mats[[elm_indx]] <- cbind(node1_input, node2_input)
          }
          
          # Do config model for within-comm
          if (i2 == i1) {
            
            degs <- nodeDegMat[nodes1, i1]
            if (sum(degs) %% 2 != 0) {
              degs[which.max(degs)] <- max(degs) + 1
            }
            edge_vec <- sample(rep.int(nodes1, times = degs))
            edgelist_mats[[elm_indx]] <- matrix(edge_vec, ncol = 2)
            
          }
          
        }
        
      }
      
    } else {
      
      edgelist_mats <- rep(list(NULL), N - 1)
      
      for (u in 1:(N - 1)) {
        
        u_nodes <- (u + 1):N
        Pi <- membership[u_nodes]
        Pj <- membership[u]
        comm_p <- P[Pi, Pj]
        pvec_u <- degrees[u_nodes] * degrees[u] / dT * comm_p
        if (max(pvec_u) > 1) pvec_u <- pvec_u / max(pvec_u)
        edges <- u_nodes[which(rbinom(N - u, 1, prob = pvec_u) == 1)]
        edgelist_mats[[u]] <- cbind(rep(u, length(edges)), edges)
        
      }
      
    }
    
    edgelist <- do.call(rbind, edgelist_mats)
    
  } else {
    
    probs <- dT * P / sum(P)
    edgelist <- generate_sbm(membership, probs, out_degs = degrees)
  }
  
  # Ordering edgelist
  unordered <- edgelist[ , 1] > edgelist[ , 2]
  edgelist[unordered, ] <- edgelist[unordered, 2:1]
  edgeid <- (edgelist[ , 1] - 1) * N + edgelist[ , 2]
  edgelist <- edgelist[order(edgeid), ]
  edgeid <- edgeid[order(edgeid)]
  edgelist <- edgelist[!duplicated(edgeid), ]
  edgelist <- edgelist[!edgelist[ , 1] == edgelist[ , 2], ]
  edgelist <- edgelist[order(edgelist[ , 1]), ]
  
  lookup <- integer(N)
  unique_nodes <- sort(unique(as.vector(edgelist)))
  lookup[unique_nodes] <- 1:length(unique_nodes)
  edgelist <- apply(edgelist, 2, function (C) lookup[C])
  
  # Making network and returning
  G <- graph.edgelist(edgelist, directed = FALSE)
  return(list(graph = G, membership = membership[unique_nodes],
              P = P, degrees = degrees, param_list = param_list))
  
}



