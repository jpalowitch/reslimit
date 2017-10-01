source("sims-code/netfuns.R")
source("full_dcsbm.R")
param_list <- make_param_list(N = 1000,
                              max_c = 100,
                              min_c = 20,
                              k = 20,
                              max_k = 50,
                              tau_1 = 2,
                              tau_2 = 1,
                              s2n = 3,
                              mu = 0.05)

set.seed(617756449)
degrees = NULL
P = NULL
membership = NULL
muversion = TRUE
type = "slow"

# Making a notes character vector
notes <- character(0)

# Extracting from paramList  
N <- ifelse(is.null(membership), param_list$N, length(membership))
max_c <- param_list$max_c
min_c <- param_list$min_c
k <- ifelse(is.null(membership), param_list$k, nrow(P))
max_k <- param_list$max_k
tau_1 <- param_list$tau_1
tau_2 <- param_list$tau_2
s2n <- param_list$s2n
mu <- param_list$mu

cat("got here 1\n")

# Making membership
if (is.null(membership)) {
  
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

cat(head(membership), "\n")

cat("got here 2\n")

if (is.null(P)) {
  # Making P matrix
  K <- max(membership)
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

cat("got here 3\n")

# Setting degrees
if (is.null(degrees)) {
  degrees <- power_law(N, tau_1, max_set = max_k, mean_set = k)
} 
dT <- sum(degrees)

# Fixing degrees step 1
if (max(degrees) > N) {
  message("adjusting degrees, since max(degrees) > N\n")
  degrees[degrees > N] <- N
}

cat("got here 4\n")
 
repeat { 
  
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
      
      cat(max_Pnum, "\n")
      
      cat("got here 4.5\n")
      
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
  cat("iteration--\n")
  cat("--", round(target_sum, 2), "vs", round(actual_sum, 2), "\n")
  if (target_sum <= actual_sum) break
}
