library(gtools)
library(gdata)
library(Matrix)
library(igraph)

# Check image.plot{fields} for a way to plot all block-plots on same page

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
  
  p_law_sample <- sample(min_set:max_set, n, replace = TRUE, prob = prob)
  return(p_law_sample)
}

# Function to make a paramList
make_param_list2 <- function (N = 5000,
                              max_c = N * 3 / 10,
                              min_c = max_c * 2 / 3,
                              k = round(sqrt(N)),
                              max_k = 3 * k,
                              tau_1 = 2,
                              tau_2 = 1,
                              mu = 0) {
  return(list("N" = N,
              "max_c" = max_c,
              "min_c" = min_c,
              "k" = k,
              "max_k" = max_k,
              "tau_1" = tau_1,
              "tau_2" = tau_2,
              "mu" = mu))
}

make_sbm <- function (param_list,
                      input_truth = NULL) {
  
  # Making a notes character vector
  notes <- character(0)
  
  # Extracting from paramList  
  N <- param_list$N
  max_c <- param_list$max_c
  min_c <- param_list$min_c
  k <- param_list$k
  max_k <- param_list$max_k
  tau_1 <- param_list$tau_1
  tau_2 <- param_list$tau_2
  mu <- param_list$mu
  
  
  if (is.null(input_truth)) {
    
    # Adding communities until amount of comms is >= o_m + 1
    # and until the sizes sum to greater than n_m
    n_m <- N
    comm_sizes <- integer(0)
    n_c <- 0
    N_comm <- sum(comm_sizes)
    while (n_c < 1 || N_comm < n_m) {
      new_draw <- power_law(1,
                            max_set = max_c,
                            min_set = min_c,
                            tau = tau_2)
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
    expanded_node_list <- unlist(lapply(1:N, 
                                        function (node) rep(node, 
                                                            n_memberships[node])))
    comm_indx_vec <- unlist(lapply(1:n_c,
                                   function (cindx) rep(cindx,
                                                        comm_sizes[cindx])))
    comm_indx_vec <- sample(comm_indx_vec)
    truth <- lapply(1:n_c, function (cindx) expanded_node_list[comm_indx_vec == cindx])
    truth <- list("communities" = truth)
    get_membership <- function (node) {
      comm_indx_vec[expanded_node_list == node]
    }
    truth$memberships <- lapply(1:N, get_membership)
    
  } else {
    truth <- input_truth
  }
  
  # Draw dVec from a power law
  dVec <- power_law(N, tau_1, max_set = max_k, mean_set = k)
  dT <- sum(dVec); pT <- sum(dVec) / N^2
  
  # Make number of edges
  nEdges <- rbinom(1, N^2, pT)
  
  # Find relative community weights
  dC <- unlist(lapply(truth$communities, function (C) sum(dVec[C])))
  pC <- dC / dT
  
  # Giving edges to communities
  cEdgeCounts <- rmultinom(1, nEdges, pC)
  
  # Find comm_vec
  nodes1 <- rep(1:N, each = N)
  nodes2 <- rep(1:N, N)
  keepers <- nodes1 < nodes2
  nodes1 <- nodes1[keepers]
  nodes2 <- nodes2[keepers]
  edgeList <- data.frame("node1" = nodes1,
                         "node2" = nodes2,
                         "weight" = rep(0, length(nodes1)))
  comm_vec <- rep(0, sum(keepers))
  for (i in 1:length(truth$communities)) {
    one_locs <- nodes1 %in% truth$communities[[i]] &
      nodes2 %in% truth$communities[[i]]
    comm_vec[one_locs] <- 1
  }
  
  n <- N + hv
  if (hv > 0) {
    bgNodes <- (N + 1):n
  } else {
    bgNodes <- integer(0)
  }
  commNodes <- 1:N
  truth$background <- bgNodes
  
  # Set phi initially
  phi <- dVec
  
  # Need target dT_C0 for community part
  # ('0' above stands for "before background is added")
  from_bg <- dVec[commNodes] * sum(dVec[bgNodes]) / dT
  target_dT_C0 <- sum(dVec[commNodes]) - sum(from_bg)
  
  
  # Calculate the pVec incorporating s_e
  pVec <- s2n_e^(comm_vec - 1) * phi[nodes1] * phi[nodes2] / dT
  dT_C_new <- sum(pVec) * 2
  
  # Calculate tol and gap
  tol <- max(pVec)
  gap <- dT_C_new / target_dT_C0
  
  # Make adjustment
  if (gap >= tol) {
    pVec <- pVec / gap
    phi[commNodes] <- phi[commNodes] / sqrt(gap)
  } else {
    pVec <- pVec / tol
    phi[commNodes] <- phi[commNodes] / sqrt(tol)
    notes <- c(notes, "gap was not above tolerance. consistency conditions probably not satisfied.")
  }
  
  # Calculate initial sVec and wVec
  sVec <- dVec^beta
  sT <- sum(sVec)
  
  # Set psi initially
  psi <- dVec^(beta - 1)
  
  # Need target sT_C0 for community part
  # ('0' above stands for "before background is added")
  from_bg_w <- sVec[commNodes] * sum(sVec[bgNodes]) / sT
  target_sT_C0 <- sum(sVec[commNodes]) - sum(from_bg_w)
  
  # Calculate wVec incorporating s_w
  wVec <- s2n_w^(comm_vec - 1) * psi[nodes1] * psi[nodes2]
  sT_C_new <- sum(wVec * pVec) * 2
  
  # Adjusting
  gap_w <- sT_C_new / target_sT_C0
  psi[commNodes] <- psi[commNodes] / sqrt(gap_w)
  wVec <- wVec / gap_w
  
  if (wtType == "exp") {
    wtSim <- function(mu, n0 = 1) 
      return(rexp(n0, 1 / mu))
  }
  
  if (wtType == "gamma") {
    wtSim <- function(mu, n0 = 1) 
      return(mu * rgamma(n0, cv^2, cv^2))
  }
  
  if (!wtType %in% c("exp", "gamma")) {
    cat("err: unrecognized wtType\n")
    return(NULL)
  }

  # If no background, need to get expectations right now
  if (hv == 0) {
    # Getting expected degrees and strengths from existing wVec, pVec
    by_nodes1 <- tapply(pVec, nodes1, sum)
    by_nodes2 <- tapply(pVec, nodes2, sum)
    Edegrees <- rep(0, N)
    to_nodes1 <- as.numeric(names(by_nodes1))
    to_nodes2 <- as.numeric(names(by_nodes2))
    Edegrees[to_nodes1] <- by_nodes1 + Edegrees[to_nodes1]
    Edegrees[to_nodes2] <- by_nodes2 + Edegrees[to_nodes2]

    by_nodes1 <- tapply(wVec * pVec, nodes1, sum)
    by_nodes2 <- tapply(wVec * pVec, nodes2, sum)
    Estrengths <- rep(0, N)
    to_nodes1 <- as.numeric(names(by_nodes1))
    to_nodes2 <- as.numeric(names(by_nodes2))
    Estrengths[to_nodes1] <- by_nodes1 + Estrengths[to_nodes1]
    Estrengths[to_nodes2] <- by_nodes2 + Estrengths[to_nodes2]
  }
  
  # Getting edges and filtering edgeList  
  edges <- rbinom(length(pVec), 1, prob = pVec)
  edgeList <- edgeList[edges == 1, ]
  pVec <- pVec[edges == 1]
  wVec <- wVec[edges == 1]
  rm(edges)
  gc()
  
  # Getting edge weights
  weights <- sapply(wVec, wtSim)
  edgeList$weight <- weights
  
  rownames(edgeList) <- NULL
  
  # Adding background, if needed
  if (hv > 0) {
    
    
    truth$background <- bgNodes

    # Getting observed in-community strengths and degrees
    adjMat <- sparseMatrix(i = edgeList$node1,
                           j = edgeList$node2,
                           x = edgeList$weight,
                           dims = c(N, N),
                           symmetric = TRUE)
    dC_O <- colSums(adjMat > 0)
    sC_O <- colSums(adjMat)
    
    # Finding adjusted target total degree
    dT_BG <- as.numeric(sum(dVec[bgNodes]))
    dT_C <- as.numeric(sum(dVec[commNodes]))
    dT_C_O <- as.numeric(sum(dC_O))
    dT_O <- sqrt((dT_BG + dT_C_O)^2 / 4 + dT_C * dT_BG) + (dT_BG + dT_C_O) / 2
    
    # Getting expected degrees for comm part
    pBG <- sum(dVec[bgNodes]) / dT_O
    EdegreesC <- dC_O + dVec[commNodes] * pBG
    
    # Finding adjusted target total strength
    sT_BG <- as.numeric(sum(sVec[bgNodes]))
    sT_C <- as.numeric(sum(sVec[commNodes]))
    sT_C_O <- as.numeric(sum(sC_O))
    sT_O <- sqrt((sT_BG + sT_C_O)^2 / 4 + sT_C * sT_BG) + (sT_BG + sT_C_O) / 2
    
    # Getting expected strengths for comm part
    wBG <- sum(sVec[bgNodes]) / sT_O
    EstrengthsC <- sC_O + sVec[commNodes] * wBG
    
    # Total expected degrees and strengths
    Estrengths <- c(EstrengthsC, sVec[bgNodes])
    Edegrees <- c(EdegreesC, dVec[bgNodes])
    
    # Making node lists
    nodes1_hv <- rep(1:n, each = hv)
    nodes2_hv <- rep(bgNodes, n)
    keepers <- nodes1_hv < nodes2_hv
    nodes1_hv <- nodes1_hv[keepers]
    nodes2_hv <- nodes2_hv[keepers]
    
    # Getting pVec_hv and wVec_hv
    pVec_hv <- Edegrees[nodes1_hv] * Edegrees[nodes2_hv] / dT_O
    wVec_hv <- Estrengths[nodes1_hv] * Estrengths[nodes2_hv] / sT_O
    
    if (max(pVec_hv) > 1) {
      notes <- c(notes, "Had to adjust pVec when adding homeless vertices")
      pVec_hv[pVec_hv > 1] <- 1
    }
    
    # Setting up hv edgelist
    edgeList_hv <- data.frame("node1" = nodes1_hv,
                              "node2" = nodes2_hv,
                              "weight" = rep(0, length(nodes1_hv)))
    
    # Simulating edges
    edges_hv <- rbinom(length(pVec_hv), 1, prob = pVec_hv)
    pVec_hv <- pVec_hv[edges_hv == 1]
    edgeList_hv <- edgeList_hv[edges_hv == 1, ]
    
    # Simulating weights
    wVec_hv <- wVec_hv[edges_hv == 1]
    weights_hv <- sapply(wVec_hv / pVec_hv, wtSim) 
    edgeList_hv$weight <- weights_hv
    
    # Concatenating edgelists
    edgeList <- rbind(edgeList, edgeList_hv)
    edgeList <- edgeList[order(edgeList$node1), ]
    
  }
  
  
  return(list("edge_list" = edgeList,
              "param_list" = param_list,
              "phi" = phi,
              "psi" = psi,
              "phi_scale" = dT,
              "Estrengths" = Estrengths,
              "Edegrees" = Edegrees,
              "truth" = truth,
              "notes" = notes))
}

sbm_summary <- function (sbm) { 
  
  edge_list <- sbm$edge_list
  param_list <- sbm$param_list
  membership <- sbm$membership
  
  N <- param_list$N
  n_c <- max(membership)
  s2n_e <- param_list$s2n
  n <- N
  commNodes <- 1:N
  
  cat("Setting up calculations...\n")
  
  # Find comm_vec
  nodes1 <- rep(1:N, each = N)
  nodes2 <- rep(1:N, N)
  keepers <- nodes1 < nodes2
  nodes1 <- nodes1[keepers]
  nodes2 <- nodes2[keepers]
  edgeList <- data.frame("node1" = nodes1,
                         "node2" = nodes2,
                         "weight" = rep(0, length(nodes1)))
  comm_vec <- rep(0, sum(keepers))
  for (i in 1:length(truth$communities)) {
    one_locs <- nodes1 %in% truth$communities[[i]] &
      nodes2 %in% truth$communities[[i]]
    comm_vec[one_locs] <- 1
  }
  
  pVec <- s2n_e^(comm_vec - 1) * phi[nodes1] * phi[nodes2] / phi_scale
  wVec <- s2n_w^(comm_vec - 1) * psi[nodes1] * psi[nodes2]
  
  if (hv > 0) {
  
    # Making node lists
    nodes1_hv <- rep(1:n, each = hv)
    nodes2_hv <- rep(bgNodes, n)
    keepers <- nodes1_hv < nodes2_hv
    nodes1_hv <- nodes1_hv[keepers]
    nodes2_hv <- nodes2_hv[keepers]
    
    # Getting pVec_hv and wVec_hv
    pVec_hv <- Edegrees[nodes1_hv] * Edegrees[nodes2_hv] / sum(Edegrees)
    wVec_hv <- Estrengths[nodes1_hv] * Estrengths[nodes2_hv] / sum(Estrengths)
    
    pVec <- c(pVec, pVec_hv)
    wVec <- c(wVec, wVec_hv / pVec_hv)
    nodes1 <- c(nodes1, nodes1_hv)
    nodes2 <- c(nodes2, nodes2_hv)
    
  }
  
  cat("Calculating expected and observed degrees/strengths...\n")
  
  # Creating extended weight expected value
  wVec_e <- wVec * pVec
  
  # Getting observed strengths and degrees

  adjMat <- sparseMatrix(i = edge_list$node1,
                         j = edge_list$node2,
                         x = edge_list$weight,
                         dims = c(N + hv, N + hv),
                         symmetric = TRUE)
  strengths <- colSums(adjMat)
  degrees <- colSums(adjMat > 0)
  
  # Calculating expected block averages
  hv_here <- as.numeric(hv > 0)
  E_edgeAvgs   <- matrix(0, n_c + hv_here, n_c + hv_here)
  E_weightAvgs <- matrix(0, n_c + hv_here, n_c + hv_here)
  edgeAvgs <- E_edgeAvgs
  weightAvgs <- E_edgeAvgs
  
  cat("Calculating expected and observed block averages...\n")
  
  for (i in 1:(n_c + hv_here)) {
    
    for (j in 1:(n_c + hv_here)) {
      
      if (i <= n_c && j <= n_c) {
        
        comm_i <- truth$communities[[i]]
        comm_j <- truth$communities[[j]]
        n_pairs <- length(comm_i) * length(comm_j)
        indx   <- which(nodes1 %in% comm_i & nodes2 %in% comm_j)
        E_edgeAvgs[i, j] <- 2 * sum(pVec[indx]) / n_pairs
        E_weightAvgs[i, j] <- 2 * sum(wVec_e[indx]) / n_pairs
        edgeAvgs[i, j] <- sum(adjMat[comm_i, comm_j] > 0) / n_pairs
        weightAvgs[i, j] <- sum(adjMat[comm_i, comm_j]) / n_pairs
        
      } else {
        
        if (i > n_c) {
          comm_i <- truth$background
        } else {
          comm_i <- truth$communities[[i]]
        } 
        
        if (j > n_c) {
          comm_j <- truth$background
        } else {
          comm_j <- truth$communities[[j]]
        }
        
        n_pairs <- length(comm_i) * length(comm_j)
        indx   <- which(nodes1 %in% comm_i & nodes2 %in% comm_j)
        E_edgeAvgs[i, j] <- sum(pVec[indx]) / n_pairs
        E_weightAvgs[i, j] <- sum(wVec_e[indx]) / n_pairs
        edgeAvgs[i, j] <- sum(adjMat[comm_i, comm_j] > 0) / n_pairs
        weightAvgs[i, j] <- sum(adjMat[comm_i, comm_j]) / n_pairs
        
        if (i > j) {
          E_edgeAvgs[i, j] <- E_edgeAvgs[j, i]
          E_weightAvgs[i, j] <- E_weightAvgs[j, i]
        }
        
        if (i == j) {
          E_edgeAvgs[i, j] <- E_edgeAvgs[i, j] * 2
          E_weightAvgs[i, j] <- E_weightAvgs[i, j] * 2
        }
        
      }
      
    }
    
  }
  
  if (hv > 0) {

    # Background diagnostics
    cat("Calculating background diagnostics...\n")
    
    bgDiagList <- rep(list(rep(list(NULL), 3)), length(truth$communities))
    
    for (i in 1:length(truth$communities)) {
      
      comm_i <- truth$communities[[i]]
      
      # Calculating stats from i to bg
      comm_stats <- rowSums(adjMat[comm_i, bgNodes])
      norm_factors <- strengths[comm_i] * sum(strengths[bgNodes]) / sum(strengths)
      comm_stats_norm <- comm_stats / norm_factors
      
      # Calculating stats from bg to i
      bg_stats <- rowSums(adjMat[bgNodes, comm_i])
      norm_factors <- sum(strengths[comm_i]) * strengths[bgNodes] / sum(strengths)
      bg_stats_norm <- bg_stats / norm_factors
  
      # Calculating stats from bg to bg
      bg_stats2 <- rowSums(adjMat[bgNodes, bgNodes])
      norm_factors <- sum(strengths[bgNodes]) * strengths[bgNodes] / sum(strengths)
      bg_stats_norm2 <- bg_stats2 / norm_factors
      
      bgDiagList[[i]][[1]] <- comm_stats_norm
      bgDiagList[[i]][[2]] <- bg_stats_norm
      bgDiagList[[i]][[3]] <- bg_stats_norm2
      
    }

  } else {
    bgDiagList <- NULL
  }
  
  # Calculating cross-block test statistics
  cat("Calculating cross-block test statistics...\n")
  
  #-----------------------------------------------------------------------
    # Variance calc
    sT <- sum(strengths)
    dT <- sum(degrees)
    suv <- strengths[edge_list$node1] * strengths[edge_list$node2] / sT
    duv <- degrees[edge_list$node1] * degrees[edge_list$node2] / dT
    duv[duv > 1] <- 1
    SSo <- sum((edge_list[,3] - suv/duv)^2)
    SSe <- sum((suv/duv)^2)
    theta <- SSo/SSe
    
    # p-value function  
    statFun <- function (B, nSet = V) {
      
      nodesToB <- nSet
      m <- length(B)
      
      sB_out <- strengths[B]
      sB_in  <- strengths[nodesToB]
      means <- sB_in * sum(sB_out) / sT    
      
      dB_out <- degrees[B]
      dB_in  <- degrees[nodesToB]
      oneB   <- rep(1, m)
      
      varListFun <- function (indx) {
        dVec   <- pmin(degrees[indx] * dB_out / dT, oneB)
        retVec <- (1 - dVec + theta) * dVec^(-1) * sB_out^2
        return (sum (retVec * strengths[indx]^2 / sT^2))
      }
      
      varList <- lapply(nodesToB, varListFun)
      vars <- unlist(varList)        
      
      stats <- rowSums(adjMat[nodesToB,B])
      norm_stats <- (stats - means) / sqrt(vars)
      return(norm_stats)
      
    }
  #-----------------------------------------------------------------------
    
  statList <- rep(list(rep(list(NULL), n_c + hv_here)), n_c + hv_here)
  statMat <- matrix(0, n_c + hv_here, n_c + hv_here)
  statMat2 <- matrix(0, n_c, n_c)
  for (i in 1:(n_c + hv_here)) {
    
    for (j in 1:(n_c + hv_here)) {
      
      cat("comms", i, "and", j, "\n")
      
      if (i <= n_c  && j <= n_c) {
      
        comm_i <- which(unlist(lapply(truth$memberships, function (m) i %in% m)))
        comm_j <- which(unlist(lapply(truth$memberships, function (m) j %in% m)))
        comm_ij <- union(comm_i, comm_j)
        
      } else {
        
        if (i > n_c) {
          
          comm_i <- truth$background
          comm_j <- which(unlist(lapply(truth$memberships, 
                                        function (m) j %in% m)))
        }
        if (j > n_c) {
          
          comm_j <- truth$background
          comm_i <- which(unlist(lapply(truth$memberships, 
                                        function (m) i %in% m)))
        }
        if (i == j) {
          comm_i <- truth$background; comm_j <- truth$background
        }
        
      }
      
      norm_stats <- statFun(comm_i, comm_j)
      norm_stats2 <- statFun(comm_ij, comm_j)
      
      statMat[i, j] <- mean(norm_stats)
      
      if (i <= n_c  && j <= n_c)
        statMat2[i, j] <- mean(norm_stats2)
      
    }
    
  }
  
  return(list("truth" = truth,
              "param_list" = sbm$param_list,
              "edge_list" = edge_list,
              "comm_vec" = comm_vec,
              "pVec" = pVec,
              "wVec" = wVec,
              "degrees" = degrees,
              "strengths" = strengths,
              "E_strengths" = Estrengths,
              "E_degrees" = Edegrees,
              "edgeAvgs" = edgeAvgs,
              "weightAvgs" = weightAvgs,
              "E_edgeAvgs" = E_edgeAvgs,
              "E_weightAvgs" = E_weightAvgs,
              "bgDiagList" = bgDiagList,
              "statMat" = statMat,
              "statMat2" = statMat2))
      
  
}

sbm_summary_plot <- function (obj, nodes = 1000) {
  
  
  if (names(obj)[1] == "truth") {
    summ <- obj
  } else {
    cat("Calculating summary...\n")
    summ <- sbm_summary(obj)
  }
  
  bes  <- summ$edgeAvgs
  bws  <- summ$weightAvgs
  Ebes <- summ$E_edgeAvgs
  Ebws <- summ$E_weightAvgs
  param_list <- summ$param_list
  comm_vec <- summ$comm_vec
  pVec <- summ$pVec
  wVec <- summ$wVec
  bgDiagList <- summ$bgDiagList
  hv <- summ$param_list$hv
  N <- summ$param_list$N
  statMat <- summ$statMat
  statMat2 <- summ$statMat2
  
  s2n_e <- param_list$s2n_e
  s2n_w <- param_list$s2n_w
  
  strengths <- summ$strengths
  degrees <- summ$degrees
  E_strengths <- summ$E_strengths
  E_degrees <- summ$E_degrees
  
  truth <- summ$truth
  edge_list <- summ$edge_list
  node_order <- c(unlist(truth$communities), truth$background)
  
  adjMat <- sparseMatrix(i = edge_list$node1,
                         j = edge_list$node2,
                         x = edge_list$weight,
                         dims = c(N + hv, N + hv),
                         symmetric = TRUE)
  
  nodes_to_plot <- sort(sample(N + hv, min(nodes, N + hv), replace = FALSE))
  plotMat <- adjMat[node_order[nodes_to_plot], node_order[nodes_to_plot]]
  
  # Plot1
  print(image(plotMat, colorkey = T))
  
  # Preparing edge/weight expectation density plots
  plot_ord <- c(which(comm_vec == 1))
  plot_ord <- c(plot_ord, setdiff(1:length(comm_vec), plot_ord))
  
  edges0 <- sample(which(comm_vec == 0), 10000, replace = FALSE)
  edges1 <- sample(which(comm_vec == 1), 10000, replace = FALSE)
  
  # Plot2
  if (hv > 0) {
    bgNodes <- (N + 1):(N + hv)
  } else {
    bgNodes <- integer(0)
  }

  cat ("Press [enter] for plot 2")
  line <- readline()

  par(mfrow = c(2, 3))
  
  hist(strengths, breaks = min(round(N / 10), 100),
       xlab = "Strength dist.")
  plot(E_strengths, strengths,
       pch = 1, xlab = "Expected strengths",
       ylab = "Observed Strengths",
       main = "Expected vs. observed strengths")
  abline(0, 1, col = "red")
  if (hv > 0) {
    points(E_strengths[bgNodes], strengths[bgNodes], col = "green")
  }
  plot(density(wVec[edges0] * s2n_w), 
       main = "w1_expv vs. s2n_w * w0_expv")
  lines(density(wVec[edges1]), col = "red")
  
  hist(degrees, breaks = min(round(N / 10), 100),
       xlab = "Degree dist.")
  plot(E_degrees, degrees,
       pch = 1, xlab = "Expected degrees",
       ylab = "Observed Degrees",
       main = "Expected vs. observed degrees")
  abline(0, 1, col = "red")
  if (hv > 0) {
    points(E_degrees[bgNodes], degrees[bgNodes], col = "green")
  }
  plot(density(pVec[edges0] * s2n_e),
       main = "e1_expv vs. s2n_w * e0_expv")
  lines(density(pVec[edges1]), col = "red")
  
  # Plot 3
  cat ("Press [enter] for plot 3.1")
  line <- readline()
  
  print(image(Matrix(Ebes, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Expected Edge Sums"))
  
  cat ("Press [enter] for plot 3.2")
  line <- readline()
  
  print(image(Matrix(bes, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Observed Edge Sums"))
  
  cat ("Press [enter] for plot 3.3")
  line <- readline()
  
  print(image(Matrix(Ebws, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Expected Weight Sums"))
  
  cat ("Press [enter] for plot 3.4")
  line <- readline()
  
  print(image(Matrix(bws, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Observed Weight Sums"))
  
  if (hv > 0) {
    # Plot 4
    
    for (i in 1:length(truth$communities)) {
      
      cat (paste0("Press [enter] for plot 4.", i))
      line <- readline()
      
      par(mfrow = c(1,3))
  
      hist(bgDiagList[[i]][[1]], xlab = "", 
           breaks = round(5 * log(length(bgDiagList[[i]][[1]]))),
           main = paste0("Bg diagnostics, comm", i, " to bg"))
      abline(v = 1, col = "red", lwd = 2)
  
      hist(bgDiagList[[i]][[2]], xlab = "",
           breaks = round(5 * log(length(bgDiagList[[i]][[1]]))),
           main = paste0("Bg diagnostics, bg to comm", i))
      abline(v = 1, col = "red", lwd = 2)
  
      hist(bgDiagList[[i]][[3]], xlab = "",
           breaks = round(5 * log(length(bgDiagList[[i]][[1]]))),
           main = paste0("Bg diagnostics, bg to bg"))
      abline(v = 1, col = "red", lwd = 2)
  
      
    }
  } else {
    cat("\n")
    cat("No bg nodes, so no plot 4\n")
    cat("\n")
  }
  
  cat ("Press [enter] for plot 5")
  line <- readline()
  
  print(image(Matrix(statMat, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Avg. z-stat from [col] nodes to [row] comm"))
  
  
  cat ("Press [enter] for plot 6")
  line <- readline()
  
  print(image(Matrix(statMat2, sparse = T),
              colorkey = T,
              cuts = 100,
              main = "Avg. z-stat from [col] nodes to [col+row] comm"))
  
  cat ("Press [enter] to turn plot off")
  line <- readline()
  
  dev.off()
  
  return("Finished")
  
}

make_oslom_run_script <- function (root_dir, 
                                   par_dirs = as.character(seq(5, 95, 5)),
                                   build_dir = "methodFiles/OSLOM2",
                                   nreps = 1) {
  
  # Getting par settings
  par_divs <- length(par_dirs)
  
  # Getting rep folder names
  rep_dirs <- as.character(1:nreps)
  
  # Initializing lines
  run_lines <- character(0)
  
  for (p in 1:par_divs) {
    
    cat("par num", p, "- ")
    
    # Setting directory
    curr_dir_p <- file.path(root_dir, par_dirs[p])
    
    for (rep in 1:nreps) {
      
      if (rep == 1)
        cat("rep: ")
      cat(rep)
      if (rep == nreps)
        cat("\n")
      
      # Setting directory
      curr_dir_p_rep <- file.path(curr_dir_p, rep_dirs[rep])
      
      run_line1 <- paste0("cd ", file.path(getwd(), curr_dir_p_rep))
      run_line2 <- paste0(file.path(getwd(), build_dir), "/oslom_undir",
                          " -f network.dat -w -singlet")
      
      # Making run string
      run_lines <- c(run_lines, run_line1, run_line2)
      
    }
    
  }
  
  writeLines(run_lines, file.path(root_dir, "oslom_run_script.txt"))
  
  return(run_lines)
  
}

save_result_dat <- function (comm_list, dir, tag) {
  
  # Getting the lines
  n_c_found <- length(comm_list)
  if (n_c_found > 0) {
    lines <- character(0)
    for (c in 1:n_c_found)
      lines <- c(lines, paste(sort(comm_list[[c]]), collapse = " "))
    writeLines(lines, file.path(dir, paste0(tag, "_comms.dat")))
  } else {
    writeLines("", file.path(dir, paste0(tag, "_comms.dat")))
  }
  
  return(tag)
  
}

read_mutual <- function (fn) {
  
  mutual_line <- readLines(fn)
  mutual_score <- as.numeric(strsplit(mutual_line, "\t")[[1]][2])
  return(mutual_score)
  
}

delete_interedges <- function (sbm) {
  
  sVec <- sbm$sd_params$sVec
  dVec <- sbm$sd_params$dVec
  truth <- sbm$truth
  hv <- sbm$param_list$hv
  edge_list <- sbm$edge_list
  N <- sbm$param_list$N
  hv <- length(truth$background)
  n_c <- length(sbm$truth$communities)
  r_e <- sbm$param_list$r_e
  r_w <- sbm$param_list$r_w
  all_strengths <- sbm$sd_params$all_strengths
  all_degrees <- sbm$sd_params$all_degrees
  
  # Making initial edgeList
  edgeList <- data.frame("node1" = rep(1:N, each = N),
                         "node2" = rep(1:N, N))
  
  
  together_vec <- integer(N^2)
  for (i in 1:n_c) {
    together_locs <- which(edgeList$node1 %in% truth$communities[[i]] &
                             edgeList$node2 %in% truth$communities[[i]])
    together_vec[together_locs] <- 1
  }
  edgeList$comm_pair <- together_vec
  
  # A function to help with indexing
  toIndex <- function (i, j, r = max(i))
    i + r * (j - 1)
  
  mat_pos <- toIndex(edgeList$node1, edgeList$node2, N)
  
  mat_pos_O <- toIndex(sbm$edge_list$node1, sbm$edge_list$node2, N)
  together_vec_O <- together_vec[match(mat_pos_O, mat_pos)]
  
  edge_list_new <- sbm$edge_list[together_vec_O == 1, ]
  
  sbm_new <- sbm
  sbm_new$edge_list <- edge_list_new
  
  return(sbm_new)
  
}

mean2mean_plot <- function (edge_list, N,
                            input_strengths = NULL,
                            input_degrees = NULL,
                            normalized = TRUE,
                            input_theta = NULL) {
  
  adjMat <- sparseMatrix(i = edge_list$node1,
                         j = edge_list$node2,
                         x = edge_list$weight,
                         dims = c(N, N),
                         symmetric = TRUE)
  
  if (is.null(input_strengths)) {
    strengths <- colSums(adjMat)
  } else {
    strengths <- input_strengths
  }
  
  if (is.null(input_degrees)) {
    degrees <- colSums(adjMat > 0)
  } else {
    degrees <- input_degrees
  }
  
  sVec <- strengths[edge_list$node1] * strengths[edge_list$node2] / 
    sum(strengths)
  dVec <- degrees[edge_list$node1] * degrees[edge_list$node2] / 
    sum(degrees)
  
  if (is.null(input_theta)) {
    SSo <- (edge_list$weight - sVec / dVec)^2
    SS <- sum((sVec / dVec)^2)
    theta <- SSo / SS
  } else {
    theta <- input_theta
  }
  
  
  means <- sVec / dVec
  stats <- edge_list$weight - means
  vars <- means^2 * input_theta
  
  colorRampFun <- colorRampPalette(c("blue", "red"))
  colPal <- colorRampFun(1000)
  maxS <- max(strengths)^2 / sum(strengths)
  minS <- min(strengths)^2 / sum(strengths)
  sProp <- (sVec - minS) / maxS
  sProp <- ceiling(sProp * 1000)
  colors <- colPal[sProp]
  
  if (normalized) {
    stats <- stats / vars
  }
  
  plot(sVec, stats, col = colors)
  abline(0, 0)
  
  return(list("stats" = stats,
              "means" = means,
              "vars" = vars))
}

check_overlap <- function (sbm, nodes_to_plot = sbm$param_list$N) {
  
  sVec <- sbm$sd_params$sVec
  dVec <- sbm$sd_params$dVec
  truth <- sbm$truth
  hv <- sbm$param_list$hv
  edge_list <- sbm$edge_list
  N <- sbm$param_list$N
  hv <- length(truth$background)
  n_c <- length(sbm$truth$communities)
  s2n_e <- sbm$param_list$s2n_e
  s2n_w <- sbm$param_list$s2n_w
  all_strengths <- sbm$sd_params$all_strengths
  all_degrees <- sbm$sd_params$all_degrees
  
  
  head(edge_list)
  
  # getting membership vector
  comm_vec <- rep(0, nrow(edge_list))
  for (c in 1:length(truth)) {
    comm_vec[edge_list$node1 %in% truth$communities[[c]] &
               edge_list$node2 %in% truth$communities[[c]]] <- 1
  }
  
  which_on <- which(unlist(lapply(truth$memberships, length)) > 1)
  on_sums <- matrix(0, length(which_on), 2)
  sn_sums <- matrix(0, N - length(which_on), 2)
  
  adjMat <- sparseMatrix(i = edge_list$node1,
                         j = edge_list$node2,
                         x = edge_list$weight,
                         dims = c(N, N),
                         symmetric = TRUE)
  
  
  for (i in 1:N) {
    
    if (i %% 10 == 0)
      cat(i, "\n")
    
    
    in_nodes <- unlist(lapply(truth$memberships[[i]], function (indx) truth$communities[[indx]]))
    if (i %in% which_on) {
      on_sums[match(i, which_on), 1] <- mean(adjMat[i, in_nodes]) / sVec[i]
      on_sums[match(i, which_on), 2] <- mean(adjMat[i, -in_nodes]) / sVec[i]
    } else {
      sn_sums[match(i, setdiff(1:N, which_on)), 1] <- mean(adjMat[i, in_nodes]) / sVec[i]
      sn_sums[match(i, setdiff(1:N, which_on)), 2] <- mean(adjMat[i, -in_nodes]) / sVec[i]
    }
    
  }
  cat("on sums:", fivenum(on_sums), "\n")
  cat("sn sums:", fivenum(sn_sums), "\n")
  
  # Finding plot order
  node_list <- rep(list(NULL), length(truth$communities) * 2 - 1)
  for (i in 1:length(truth$communities)) {
    ol_prev <- integer(0)
    if (i < length(truth$communities)) {
      ol_flag <- function (ms) {sum(c(i, i + 1) %in% ms) == 2}
      ol_nodes <- which(unlist(lapply(truth$memberships, ol_flag)))
      node_list[[2 * i - 1]] <- setdiff(truth$communities[[i]],
                                        union(ol_nodes, ol_prev))
      node_list[[2 * i]] <- ol_nodes
      ol_prev <- ol_nodes
    } else {
      node_list[[2 * i - 1]] <- setdiff(truth$communities[[i]],
                                        ol_prev)
    }
  }
  node_vec <- unlist(node_list)
  
  plot_node_index <- sort(sample(1:N, nodes_to_plot, replace = FALSE))
  plot_these <- node_vec[plot_node_index]
  plotMat <- adjMat[plot_these, plot_these]
  
  print(image(plotMat, colorkey = T))
  
}



