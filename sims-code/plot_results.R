Args = commandArgs(trailingOnly=TRUE)

# Compute times in log values?
log_time <- TRUE

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

to_run <- as.numeric(Args)

library(ggplot2)
library(reshape2)

for (exper in exper_names[to_run]) {
  
  # Loading/setting parameters
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  # Finding methNames
  dummy_dir <- file.path(rootdir, 1)
  fns <- list.files(dummy_dir)
  methFiles <- fns[grepl("results", fns)]
  methNames <- unique(sapply(methFiles, function (c) strsplit(c, "_")[[1]][1]))
  rm(dummy_dir, fns, methFiles)
  
  # Loading results matrices
  timers <- read.table(file.path(rootdir, "timer_mat.txt"), 
                       row.names = 1, header = FALSE)
  onmis  <- read.table(file.path(rootdir, "onmi_mat.txt"), 
                       row.names = 1, header = FALSE)
  CIBs <- read.table(file.path(rootdir, "CIB_mat.txt"), 
                       row.names = 1, header = FALSE)
  
  # Formatting plot data
  timer_plotdf <- melt(t(timers))
  onmis_plotdf <- melt(t(onmis))
  CIBs_plotdf <- melt(t(CIBs))
  levels(timer_plotdf$Var1) <- levels(onmis_plotdf$Var1) <- 
    levels(CIBs_plotdf$Var1) <- ps
  timer_plotdf$Var1 <- as.numeric(levels(timer_plotdf$Var1))
  onmis_plotdf$Var1 <- as.numeric(levels(onmis_plotdf$Var1))
  CIBs_plotdf$Var1 <- as.numeric(levels(CIBs_plotdf$Var1))
  
  # Plotting and saving
  timer_str <- ifelse(log_time, "log10(seconds + 1)", "seconds")
  tp <- ggplot(timer_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab(timer_str) + 
    ggtitle(paste0(exper, " Timing")) + 
    guides(colour = guide_legend(title = "Method"))
  op <- ggplot(onmis_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("ONMI") + 
    ggtitle(paste0(exper, " ONMI")) + 
    guides(colour = guide_legend(title = "Method"))
  cp <- ggplot(CIBs_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("CIB") + 
    ggtitle(paste0(exper, " CIB")) + 
    guides(colour = guide_legend(title = "Method"))
  ggsave(file.path(rootdir, "times.png"), tp)
  ggsave(file.path(rootdir, "onmis.png"), op)
  ggsave(file.path(rootdir, "CIBs.png"), cp)
  
  # Making pairwise ONMI plots
  non_rmods <- setdiff(levels(timer_plotdf$Var2), "rmodbLouv")
  for (method in non_rmods) {
    meth_plotdf <- subset(onmis_plotdf, onmis_plotdf$Var2 %in% c(method, "rmodbLouv"))
    op <- ggplot(meth_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
          geom_line() + xlab(pname) + ylab("ONMI") + 
          ggtitle(paste0(exper, " ONMI")) + ylim(c(0, 1)) +
          guides(colour = guide_legend(title = "Method"))
    ggsave(file.path(rootdir, paste0("rmod_vs_", method, ".png")))
  }
  
}
