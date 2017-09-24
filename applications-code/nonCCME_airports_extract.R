source("slpaRead.R")
source("osRead.R")

# Making month char list
months = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

##### Reads everything from results files and plots

source("applications-code/plotPortComms.R")
saveDir = "applications-results/airports"
monthNames = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
fnList = c(monthNames,"year")

### Specify whether or not to include plot titles
plotTitles = FALSE

### Set filter size (i.e. minimum cluster size)
filtSize = 1

## Specify years
load("applications-results/airports/data/dataInfo.RData")
yL <- lastYear - firstYear + 1


#######################################
#######################################

for (y in 1:yL) {
  
  yearNum <- firstYear + y - 1
  curr_dir <- file.path(saveDir,yearNum)
  allFiles <- list.files(curr_dir, full.names = TRUE)
  outputFiles <- allFiles[grepl("output",allFiles)]
  outputFiles <- outputFiles[grepl(".RData",outputFiles)]
  
  icpmFiles <- allFiles[grepl("icpm", allFiles)]  
  
  for (fn in icpmFiles) {
    
    # Finding which month
    month <- strsplit(fn, "_")[[1]][3]
    
    # Getting SLPA results
    results <- slpaRead(fn)
    save(results, file = file.path(curr_dir,
                                   paste0("output_slpa_",
                                          month,
                                          ".RData")))
    
    # Getting OSLOM results
    dat_fn <- paste0(yearNum, "_", month, ".dat")
    results <- osRead(file.path("applications-results/airports/OSLOM2",
                                paste0(dat_fn, "_oslo_files"),
                                "tp"))
    save(results, file = file.path(curr_dir,
                                   paste0("output_oslom_",
                                          month,
                                          ".RData")))
    
    
    
  }
  
}

