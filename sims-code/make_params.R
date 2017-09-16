sim_res_dir <- "sims-results"

################################################################################

if (!dir.exists(sim_res_dir))
  dir.create(sim_res_dir)

exper_names <- NULL

################################################################################

# twocliq with increasing nS
exper_name <- "twocliq_nC=8_increase_nS"
ps <- seq(100, 1000, 100)
pname <- "nS"
nsims <- 1
extra_pars <- list(nC = 8)
make_code <- "Gp <<- twocliq(ps[p], extra_pars$nC, 0.5)"
make_type <- "R_igraph"
truth_code <- "comms <<- list(1:ps[p], (ps[p] + 1):(ps[p] + extra_pars$nC),
                    (ps[p] + extra_pars$nC + 1):(ps[p] + 2 * extra_pars$nC))"
truth_type <- "Manual"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# LFR with increasing mu
exper_name <- "LFR_N=1000_B_20-100_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- list(N = "-N 1000",
                   k = "-k 20", maxk = "-maxk 50",
                   minc = "-minc 20", maxc = "-maxc 100")
make_code <- c("system2('./../../../binary_networks/benchmark',
                        paste(extra_pars$N, extra_pars$k, extra_pars$maxk,
                              '-mu', ps[p], extra_pars$minc, extra_pars$maxc))",
               "file.copy('network.dat', paste0(i, '.dat'))")
make_type <- "System"
truth_type <- "LFR"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# LFR with increasing mu, less harsh degrees
exper_name <- "LFR_N=1000_B_20-100_D_150-200_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- list(N = "-N 1000",
                   k = "-k 150", maxk = "-maxk 200",
                   minc = "-minc 20", maxc = "-maxc 100")
make_code <- c("system2('./../../../binary_networks/benchmark',
                        paste(extra_pars$N, extra_pars$k, extra_pars$maxk,
                              '-mu', ps[p], extra_pars$minc, extra_pars$maxc))",
               "file.copy('network.dat', paste0(i, '.dat'))")
make_type <- "System"
truth_type <- "LFR"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# SBM with increasing mu
exper_name <- "SBM_N=1000_B_20-100_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- NULL
make_code <- c("param_list <<- make_param_list(N = 1000, max_c = 100, min_c = 20, k = 20, max_k = 50, mu = ps[p])",
               "Gobj <<- DCSBM(param_list)",
               "Gp <<- Gobj$graph")
truth_code <- c("K <<- max(Gobj$membership)",
                "comms <<- lapply(1:K, function (i) which(Gobj$membership == i))")
make_type <- "R_igraph"
truth_type <- "Manual"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# SBM with increasing mu, less harsh degrees
exper_name <- "SBM_N=1000_B_20-100_D_150-200_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- NULL
make_code <- c("param_list <<- make_param_list(N = 1000, max_c = 100, min_c = 20, k = 150, max_k = 200, mu = ps[p])",
               "Gobj <<- DCSBM(param_list)",
               "Gp <<- Gobj$graph")
truth_code <- c("K <<- max(Gobj$membership)",
                "comms <<- lapply(1:K, function (i) which(Gobj$membership == i))")
make_type <- "R_igraph"
truth_type <- "Manual"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# SBM0 with increasing mu
exper_name <- "SBM0_N=1000_B_20-100_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- NULL
make_code <- c("param_list <<- make_param_list(N = 1000, max_c = 100, min_c = 20, k = 20, max_k = 50, mu = ps[p])",
               "Gobj <<- DCSBM(param_list, type = 'slow')",
               "Gp <<- Gobj$graph")
truth_code <- c("K <<- max(Gobj$membership)",
                "comms <<- lapply(1:K, function (i) which(Gobj$membership == i))")
make_type <- "R_igraph"
truth_type <- "Manual"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

# SBM0 with increasing mu, less harsh degrees
exper_name <- "SBM0_N=1000_B_20-100_D_150-200_increase_mu"
ps <- seq(0.05, 0.95, 0.05)
pname <- "mu"
nsims <- 3
extra_pars <- NULL
make_code <- c("param_list <<- make_param_list(N = 1000, max_c = 100, min_c = 20, k = 150, max_k = 200, mu = ps[p])",
               "Gobj <<- DCSBM(param_list, type = 'slow')",
               "Gp <<- Gobj$graph")
truth_code <- c("K <<- max(Gobj$membership)",
                "comms <<- lapply(1:K, function (i) which(Gobj$membership == i))")
make_type <- "R_igraph"
truth_type <- "Manual"
save(ps, pname, nsims, extra_pars, make_code, make_type, truth_code, truth_type,
     file = file.path(sim_res_dir, paste0("pars_", exper_name, ".RData")))

exper_names <- c(exper_names, exper_name)

################################################################################

writeLines(exper_names, con = file.path(sim_res_dir, "exper_names.txt"))
writeLines(sim_res_dir, con = "sim_res_dir.txt")
