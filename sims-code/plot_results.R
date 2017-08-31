library(ggplot2)
library(reshape2)

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))

for (exper in exper_names[-1]) {
  
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
  
  # Formatting plot data
  timer_plotdf <- melt(t(timers))
  onmis_plotdf <- melt(t(onmis))
  levels(timer_plotdf$Var1) <- levels(onmis_plotdf$Var1) <- ps
  timer_plotdf$Var1 <- as.numeric(levels(timer_plotdf$Var1))
  onmis_plotdf$Var1 <- as.numeric(levels(onmis_plotdf$Var1))
  
  # Plotting and saving
  tp <- ggplot(timer_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("Time (sec)") + 
    ggtitle(paste0(exper, " Timing")) + 
    guides(colour = guide_legend(title = "Method"))
  op <- ggplot(onmis_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("ONMI") + 
    ggtitle(paste0(exper, " ONMI")) + 
    guides(colour = guide_legend(title = "Method"))
  ggsave(file.path(rootdir, "times.png"), tp)
  ggsave(file.path(rootdir, "onmis.png"), op)
  
  
}
