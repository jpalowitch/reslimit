source("netfuns.R")

# Remake sims or just write scripts?
remake_sims <- TRUE

if (!dir.exists("sims2"))
  dir.create("sims2")

set.seed(12345)

# er with increasing n
rootdir <- "sims2/er_increase_n"
if (!dir.exists(rootdir))
  dir.create(rootdir)
ps <- seq(100, 500, 100)
pname <- "n"
nsims <- 1
cscor_lines <- NULL

for (p in 1:length(ps)) {
  
  # Setting directory
  curr_dir <- file.path(rootdir, p)
  if (!dir.exists(curr_dir))
    dir.create(curr_dir)
  
  for (i in 1:nsims) {
    
    # Setting seed
    seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
    if (!file.exists(seedfile)) {
      seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
      writeLines(seedpi, con = seedfile)
    }
    
    # Making and saving data
    set.seed(as.numeric(readLines(seedfile)))
    Gp <- erdos.renyi.game(n = ps[p], p = 0.5)
    dat_fn <- paste0(i, ".dat")
    write.table(get.edgelist(Gp), sep = "\t",
                file = file.path(curr_dir, dat_fn),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    save(Gp, file = file.path(curr_dir, paste0(i, ".RData")))
    
    # Writing cb-signi lines
    writeLines(paste(c(1:50), collapse = " "), con = file.path(curr_dir, "comm1.dat"))
    writeLines(paste(c(1:(ps[p] / 2)), collapse = " "), con = file.path(curr_dir, "comm2.dat"))
    cscor_add1 <- paste("{ time",
                        "cb-signi/compare", "-f",  file.path(curr_dir, dat_fn),
                        "-c", file.path(curr_dir, "comm1.dat"),
                        "-t 0.01", "-nobcore",
                        "; } 2>", file.path(curr_dir, "comm1_time.txt"))
    cscor_add2 <- paste("{ time",
                        "cb-signi/compare", "-f",  file.path(curr_dir, dat_fn),
                        "-c", file.path(curr_dir, "comm2.dat"),
                        "-t 0.01", "-nobcore",
                        "; } 2>", file.path(curr_dir, "comm2_time.txt"))
    cscor_lines <- c(cscor_lines, cscor_add1, cscor_add2)
    
  }
  
}

writeLines(cscor_lines, con = "cscor_lines_n.txt")


# er with increasing p
rootdir <- "sims2/er_increase_p"
if (!dir.exists(rootdir))
  dir.create(rootdir)
ps <- seq(0.1, 0.5, 0.1)
pname <- "p"
nsims <- 1
cscor_lines <- NULL

for (p in 1:length(ps)) {
  
  # Setting directory
  curr_dir <- file.path(rootdir, p)
  if (!dir.exists(curr_dir))
    dir.create(curr_dir)
  
  for (i in 1:nsims) {
    
    # Setting seed
    seedfile <- file.path(curr_dir, paste0(i, "_seed.txt"))
    if (!file.exists(seedfile)) {
      seedpi <- paste(sample(0:9, 9, replace = TRUE), collapse = "")
      writeLines(seedpi, con = seedfile)
    }
    
    # Making and saving data
    set.seed(as.numeric(readLines(seedfile)))
    Gp <- erdos.renyi.game(n = 500, p = ps[p])
    dat_fn <- paste0(i, ".dat")
    write.table(get.edgelist(Gp), sep = "\t",
                file = file.path(curr_dir, dat_fn),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    save(Gp, file = file.path(curr_dir, paste0(i, ".RData")))
    
    # Writing cb-signi lines
    writeLines(paste(c(1:50), collapse = " "), con = file.path(curr_dir, "comm1.dat"))
    writeLines(paste(c(1:(50 * p)), collapse = " "), con = file.path(curr_dir, "comm2.dat"))
    cscor_add1 <- paste("{ time",
                        "cb-signi/compare", "-f",  file.path(curr_dir, dat_fn),
                        "-c", file.path(curr_dir, "comm1.dat"),
                        "-t 0.01", "-nobcore",
                        "; } 2>", file.path(curr_dir, "comm1_time.txt"))
    cscor_add2 <- paste("{ time",
                        "cb-signi/compare", "-f",  file.path(curr_dir, dat_fn),
                        "-c", file.path(curr_dir, "comm2.dat"),
                        "-t 0.01", "-nobcore",
                        "; } 2>", file.path(curr_dir, "comm2_time.txt"))
    cscor_lines <- c(cscor_lines, cscor_add1, cscor_add2)
    
  }
  
}

writeLines(cscor_lines, con = "cscor_lines_p.txt")
