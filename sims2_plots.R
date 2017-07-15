dir.create("sims2/plots")

rootdir <- "sims2/er_increase_n"
npar <- length(list.files(rootdir))
fixedCtimers <- numeric(npar)
grownCtimers <- numeric(npar)

for (i in 1:npar) {
  
  curr_dir <- file.path(rootdir, i)
  t1 <- as.character(read.table(file.path(curr_dir, "comm1_time.txt"))$V2[1])
  t2 <- as.character(read.table(file.path(curr_dir, "comm2_time.txt"))$V2[1])
  ms1 <- strsplit(t1, "m")[[1]]
  ms2 <- strsplit(t2, "m")[[1]]
  ms1[2] <- substr(ms1[2], 1, 5)
  ms2[2] <- substr(ms2[2], 1, 5)
  fixedCtimers[i] <- as.numeric(ms1[1]) * 60 + as.numeric(ms1[2])
  grownCtimers[i] <- as.numeric(ms2[1]) * 60 + as.numeric(ms2[2])
  
}

png("sims2/plots/er_increase_n1.png")
plot(seq(100, 500, 100), fixedCtimers, ylab = "sec", xlab = "n",
     main = "ER, p = 0.5, |C| = 50")
dev.off()

png("sims2/plots/er_increase_n2.png")
plot(seq(100, 500, 100), grownCtimers, ylab = "sec", xlab = "n",
     main = "ER, p = 0.5, |C| = n / 2")
dev.off()

rootdir <- "sims2/er_increase_p"
npar <- length(list.files(rootdir))
fixedCtimers <- numeric(npar)
grownCtimers <- numeric(npar)

for (i in 1:npar) {
  
  curr_dir <- file.path(rootdir, i)
  t1 <- as.character(read.table(file.path(curr_dir, "comm1_time.txt"))$V2[1])
  t2 <- as.character(read.table(file.path(curr_dir, "comm2_time.txt"))$V2[1])
  ms1 <- strsplit(t1, "m")[[1]]
  ms2 <- strsplit(t2, "m")[[1]]
  ms1[2] <- substr(ms1[2], 1, 5)
  ms2[2] <- substr(ms2[2], 1, 5)
  fixedCtimers[i] <- as.numeric(ms1[1]) * 60 + as.numeric(ms1[2])
  grownCtimers[i] <- as.numeric(ms2[1]) * 60 + as.numeric(ms2[2])
  
}

png("sims2/plots/er_increase_p1.png")
plot(seq(0.1, 0.5, 0.1), fixedCtimers, ylab = "sec", xlab = "p",
     main = "ER, n = 500, |C| = 50")
dev.off()

png("sims2/plots/er_increase_p2.png")
plot(seq(0.1, 0.5, 0.1), grownCtimers, ylab = "sec", xlab = "p",
     main = "ER, n = 500, |C| = 250")
dev.off()
  
  
  
  