mod_sig <- function (fn = NULL, edgelist = NULL, mod_type = "louvain", nsims = 1000) {
  
  # Loading network file and creating graph object
  if (is.null(edgelist)) {
    edgelist <- as.matrix(read.table(fn))
  }
  G <- graph.edgelist(as.matrix(edgelist), directed = FALSE)
  n <- length(V(G))
  
  # Running optimization
  if (mod_type == "louvain") {
    res <- cluster_louvain(G)
    comms <- lapply(1:max(res$membership), 
                  function (i) which(res$membership == i))
  } else {
    stop("mod_type unsupported\n")
  }
  original_mod <- tail(res$modularity, 1)
  degs <- degree(G)
  
  # Running configuration model & mod max many times
  disp_count <- 0
  mods <- numeric(nsims)
  for (i in 1:nsims) {
    
    if (i > ceiling(nsims / 10 * disp_count)) {
      cat(paste0(disp_count * 10, '%'), "-- ")
      disp_count <- disp_count + 1
    }
  
    cm_G <- config_model(degs, no_multiple = TRUE)
    if (mod_type == "louvain") {
      cm_res <- cluster_louvain(cm_G)
    }
    mods[i] <- tail(cm_res$modularity, 1)
  }
  
  return(list(original_mod = original_mod,
              mods = mods,
              pval = 1 - sum(mods <= original_mod) / nsims))
}