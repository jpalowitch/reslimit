plot_hierarchies <- function (resolution_obj, fn) {
  
  suppressMessages(library(ggplot2))
  library(reshape)
  library(grid)
  library(gridExtra)
  library(scales)
  
  commsizes <- resolution_obj$Ks
  #commsizes <- commsizes
  consecvis <- resolution_obj$VIs
  
  p1_df <- data.frame(Resolution_Parameter = rep(resolution_obj$res_pars, 2),
                      ydata = c(commsizes, consecvis),
                      Stat = c(rep("Community_Count", length(resolution_obj$res_pars)),
                               rep("VI", length(resolution_obj$res_pars))))
  VI_range <- range(consecvis)
  CS_range <- range(commsizes)
  CS_norm <- (commsizes - CS_range[1]) / (CS_range[2] - CS_range[1])
  CS_norm <- CS_norm * (VI_range[2] - VI_range[1])
  p1_df$ydata[p1_df$Stat == "Community_Count"] <- CS_norm
  axis_map <- function (x) {
    x / (VI_range[2] - VI_range[1]) * (CS_range[2] - CS_range[1]) + CS_range[1]
  }
  
  
  p1 <- ggplot(p1_df, aes(x = Resolution_Parameter, y = ydata, colour = Stat)) + 
    geom_line() + scale_y_continuous(sec.axis = sec_axis(~axis_map(.), "log10(# of communities)",
                                                         breaks = pretty_breaks(n = 6))) + 
    ylab("Consecutive_VI") + theme(legend.position = c(0.8, 0.2)) + 
    ggtitle("RBConfig result statistics")
  
  
  p1 <- ggplot(p1_df, aes(x = Resolution_Parameter, y = ydata, colour = Stat)) + 
    geom_line() + scale_y_continuous(trans = 'log10',
                                     breaks = trans_breaks('log10', function(x) 10^x),
                                     labels = trans_format('log10', math_format(10^.x)),
                        sec.axis = sec_axis(~axis_map(.), "log10(# of communities)",
                                                         breaks = pretty_breaks(n = 6))) + 
    ylab("Consecutive_VI") + theme(legend.position = c(0.8, 0.2)) + 
    ggtitle("RBConfig result statistics")

    
  p2_df <- melt(t(resolution_obj$part_mat))
  p2_df$X1 <- resolution_obj$res_pars[p2_df$X1]
  p2_df$X2 <- factor(p2_df$X2)
  p2_df$value <- -log10(p2_df$value)
  names(p2_df) <- c("Resolution_Parameter", "Partition_Tier", "VI")
  p2 <- ggplot(p2_df, aes(x = Resolution_Parameter, 
                          y = VI, 
                          colour = Partition_Tier)) + 
    geom_line() + theme(legend.position = c(0.8, 0.8)) + 
    ylab("-log10(VI to RBConfig result)") + 
    scale_y_continuous(sec.axis = sec_axis(~identity(.), "",
                                                         breaks = pretty_breaks(n = 6))) + 
    ggtitle("Correspondence of recursive-mod results to RBConfig")
  
  ggsave(fn, grid.arrange(p1, p2, nrow = 2),
         height = 10, width = 6)
  
}
                    