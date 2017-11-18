Args = commandArgs(trailingOnly=TRUE)

# Compute times in log values?
log_time <- TRUE

sim_res_dir <- readLines("sim_res_dir.txt")
exper_names <- readLines(file.path(sim_res_dir, "exper_names.txt"))
to_run <- as.numeric(Args)

library(ggplot2)
library(reshape2)
library(RColorBrewer)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

methNames <- c("rmodbLouv", "infomap", "louvain", "modbpy", "oslom", "RBres", "sbm", "essc")
plot_names <- c("PMrecursion", "Infomap", "Louvain", "ZMrecursion", "OSLOM", "RBresolution", "GraphToolSBM", "ESSC")
colPal <- gg_color_hue(length(methNames))
my_theme <- theme(axis.text.x = element_text(size = 15),
                  axis.title.x = element_text(size = 15),
                  axis.text.y = element_text(size = 15),
                  axis.title.y = element_text(size = 15),
                  title = element_text(size = 20),
                  legend.text = element_text(size = 13),
                  legend.title = element_text(size = 15))
pw = 8
ph = 7

for (exper in exper_names[to_run]) {
  
  cat("doing exper", exper, "which is number", match(exper, exper_names), "\n")

  
  # Loading/setting parameters
  rootdir <- file.path(sim_res_dir, exper)
  load(file.path(sim_res_dir, paste0("pars_", exper, ".RData")))
  
  # Finding methNames
  #dummy_dir <- file.path(rootdir, 1)
  #fns <- list.files(dummy_dir)
  #methFiles <- fns[grepl("results", fns)]
  #methNames <- unique(sapply(methFiles, function (c) strsplit(c, "_")[[1]][1]))
  #rm(dummy_dir, fns, methFiles)
  
  # Loading results matrices
  timers <- read.table(file.path(rootdir, "timer_mat.txt"), 
                       row.names = 1, header = FALSE)
  onmis  <- read.table(file.path(rootdir, "onmi_mat.txt"), 
                       row.names = 1, header = FALSE)
  CIBs <- read.table(file.path(rootdir, "CIB_mat.txt"), 
                       row.names = 1, header = FALSE)
  
  # Formatting plot data
  timer_plotdf <- melt(t(timers[methNames, ]))
  onmis_plotdf <- melt(t(onmis[methNames, ]))
  CIBs_plotdf <- melt(t(CIBs[methNames, ]))
  levels(timer_plotdf$Var1) <- levels(onmis_plotdf$Var1) <- 
    levels(CIBs_plotdf$Var1) <- ps
  levels(timer_plotdf$Var2) <- levels(onmis_plotdf$Var2) <- 
    levels(CIBs_plotdf$Var2) <- plot_names
  timer_plotdf$Var1 <- as.numeric(levels(timer_plotdf$Var1))
  onmis_plotdf$Var1 <- as.numeric(levels(onmis_plotdf$Var1))
  CIBs_plotdf$Var1 <- as.numeric(levels(CIBs_plotdf$Var1))
  
  # Plotting and saving
  timer_str <- ifelse(log_time, "log10(seconds + 1)", "seconds")
  tp <- ggplot(timer_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab(timer_str) + 
    #ggtitle(paste0(exper, " Timers")) + ylim(c(0, 1)) +
    ggtitle("Timing Results") +
    guides(colour = guide_legend(title = "Method")) + my_theme
  op <- ggplot(onmis_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("NMI") + 
    #ggtitle(paste0(exper, " CIB")) + ylim(c(0, 1)) +
    ggtitle("Performance Results (NMI)") + ylim(c(0, 1)) +
    guides(colour = guide_legend(title = "Method")) + my_theme
  cp <- ggplot(CIBs_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
    geom_line() + xlab(pname) + ylab("CIB") + 
    #ggtitle(paste0(exper, " ONMI")) + ylim(c(0, 1)) +
    ggtitle("Performance Results (CIB)") + ylim(c(0, 1)) +
    guides(colour = guide_legend(title = "Method")) + my_theme
  ggsave(file.path(rootdir, "times.png"), tp, width = pw + 2, height = ph)
  ggsave(file.path(rootdir, "onmis.png"), op, width = pw + 2, height = ph)
  ggsave(file.path(rootdir, "CIBs.png"), cp, width = pw + 2, height = ph)
  
  # Making pairwise ONMI plots
  non_rmods <- setdiff(levels(timer_plotdf$Var2), "PMrecursion")
  for (method in non_rmods) {
    meth_plotdf <- subset(onmis_plotdf, onmis_plotdf$Var2 %in% c(method, "PMrecursion"))
    col_ord <- match(unique(meth_plotdf$Var2), plot_names)
    cols = colPal[col_ord]
    op <- ggplot(meth_plotdf, aes(x = Var1, y = value, colour = Var2)) + 
          geom_line() + xlab(pname) + ylab("NMI") + 
          #ggtitle(paste0(exper, " ONMI")) + ylim(c(0, 1)) +
          ggtitle(paste0("PMrecursion vs. ", method)) + ylim(c(0, 1)) +
          guides(colour = guide_legend(title = "Method")) + 
          my_theme + scale_colour_manual(values = cols)
    ggsave(file.path(rootdir, paste0("rmod_vs_", methNames[plot_names == method], ".png")),
           width = pw, height = ph)
  }
  
}
