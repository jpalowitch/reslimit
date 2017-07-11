source("netfuns.R")

if (!dir.exists("sims"))
  dir.create("sims")

set.seed(12345)

# twocliq with increasing nS
rootdir <- "sims/twocliq_nS"
if (!dir.exists(rootdir))
  dir.create(rootdir)
ps <- round(50 * 3^(seq(0, 4, 0.5)))
pname <- "nS"
nsims <- 3
nC <- 8
oslom_lines <- NULL

for (p in 1:length(ps)) {
  
  # Setting directory
  curr_dir <- file.path(rootdir, p)
  if (!dir.exists(curr_dir))
    dir.create(curr_dir)
  oslom_lines <- c(oslom_lines, paste("cd", curr_dir))
  
  for (i in 1:nsims) {
    
    # Setting seed
    seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
    if (!file.exists(seedfile)) {
      seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
      writeLines(seedpi, con = seedfile)
    }
    
    # Making and saving data
    set.seed(as.numeric(readLines(seedfile)))
    Gp <- twocliq(ps[p], nC, 0.5)
    dat_fn <- paste0(i, ".dat")
    write.table(get.edgelist(Gp), 
                file = file.path(curr_dir, dat_fn),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    save(Gp, file = file.path(curr_dir, paste0(i, ".RData")))
    
    # Making truth
    comms <- list(1:ps[p], (ps[p] + 1):(ps[p] + nC),
                  (ps[p] + nC + 1):(ps[p] + 2 * nC))
    truth_strings <- unlist(lapply(comms, paste, collapse = " "))
    writeLines(truth_strings,
               con = file.path(curr_dir, paste0(i, "_truth.dat")))
    
    # Writing OSLOM lines
    oslom_add <- paste("./../../../OSLOM2/oslom_undir -uw -f", 
                       dat_fn, "-singlet", "-fast")
    oslom_lines <- c(oslom_lines, oslom_add)
    
  }
  
  oslom_lines <- c(oslom_lines, paste("cd", "../../../"))
  
}

save(ps, nsims, pname, file = file.path(rootdir, "pars.RData"))

writeLines(oslom_lines, con = file.path(rootdir, "oslom_lines.txt"))
