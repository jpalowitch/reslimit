#	Setup ------------------------------------------------------------------------

library(gdata)
airDir <- "applications-results/airports/data"
webRoot <- "http://stats.johnpalowitch.com/data-sets/CCME_analyses/airports"
firstYear <- 2009
lastYear <- 2015
readInData <- TRUE

# ------------------------------------------------------------------------------

web_files <- readLines("applications-code/airport-https.txt")
if (length(web_files) != lastYear - firstYear + 1)
  stop('year span, web_file length mismatch\n')
web_addresses <- file.path(webRoot, firstYear:lastYear, web_files)
names(web_addresses) <- as.character(firstYear:lastYear)
if (!dir.exists(airDir))
  dir.create(airDir, recursive = TRUE)

years <- seq(firstYear,lastYear,1)
yL <- length(years)

if (readInData) {
  
  dataList <- rep(list(NULL), yL)
  
  for (i in 1:yL) {
  
  	cat(i, "\n")
  
  	## Getting raw table
  	yr <- years[i]
  	iFold <- file.path(airDir, yr)
  	fn <- list.files(iFold)
  	datTab <- read.csv(web_addresses[as.character(yr)],
  	                   header = TRUE,
  	                   stringsAsFactors = FALSE)
  
  	cat("months\n")	
  
  	## Creating month list
  	maxMonth <- max(datTab$MONTH)
  	monthList <- rep(list(NULL), maxMonth)
  	for (j in 1:maxMonth) {
  
  		cat(j, "\n")
  
  		monthRows <- datTab$MONTH == j
  		monthList[[j]]$dest <- datTab$DEST[monthRows]
  		monthList[[j]]$org <- datTab$ORIGIN[monthRows]
  		monthList[[j]]$count <- datTab$PASSENGERS[monthRows]
  		
  	}
  
  	dataList[[i]]$months <- monthList
  	dataList[[i]]$dest <- datTab$DEST
  	dataList[[i]]$org <- datTab$ORIGIN
  
  }

  save(dataList, file = file.path(airDir, "rawData.RData"))

} else {
  load(file.path(airDir, "rawData.RData"))
}

# Getting master code list, making skeleton adjMat
allCodes <- character(0)
for (i in 1:yL) {
	allCodes <- union(allCodes,
	                  union(unique(dataList[[i]]$dest),
	                        unique(dataList[[i]]$org)
	                        )
	                  )
}
allCodes <- sort(allCodes)
N <- length(allCodes)
adjMat0 <- matrix(0, N, N)
rownames(adjMat0) <- allCodes
colnames(adjMat0) <- allCodes

# Making month adjLists
adjList <- rep(list(NULL),yL)
names(adjList) <- sapply(seq(firstYear, lastYear, 1), toString)
for (i in 1:yL) {

	cat("#------------\n", i, "\n#------------\n")

	maxMonth <- length(dataList[[i]]$months)
	monthList <- rep(list(NULL), maxMonth)
	for (j in 1:maxMonth) {

		cat("--#----------\n", j, "\n--#----------\n")

		adjMat <- adjMat0
		M <- length(dataList[[i]]$months[[j]]$count)
		for (entry in 1:M) {
		  
			if (entry%%1 == 0)
				cat("entry ", entry, "\n")
		  
		  enter_this <- dataList[[i]]$months[[j]]$count[entry]
		  
			adjMat[dataList[[i]]$months[[j]]$org[entry],
			       dataList[[i]]$months[[j]]$dest[entry]] <- enter_this
			
		}

		lowerTriangle(adjMat) <- lowerTriangle(t(adjMat))
		diag(adjMat) <- 0
	
		monthList[[j]] <- adjMat
		rm(adjMat)
	}

	adjList[[i]]$months <- monthList
}

save(adjList, file = file.path(airDir, "adjData.RData"))

# Analyzing by month,year
# Want to look at:
#	Number of zero degrees
#	Max degree
#	Max strength
#	Overall strength
#	Sparsity

colNames = c("nZero","mDeg","mS","st","spars")

## Getting month count and making stat mats
monthCount = 0
for(i in 1:yL)
	monthCount = monthCount + length(adjList[[i]]$months)
monthMat = matrix(0,monthCount,5)
yearMat = matrix(0,yL,5)
colnames(monthMat) = colNames
colnames(yearMat) = colNames


for(i in 1:yL){

	cat("#####\n",i,"\n#####\n")

	yrAdj = adjMat0

	maxMonth = length(adjList[[i]]$months)
	for(j in 1:maxMonth){
		monthAdj = adjList[[i]]$months[[j]]
		yrAdj = yrAdj + monthAdj

		## Calc stats
		degrees = colSums(monthAdj>0)
		nZeroDegs = sum(degrees==0)
		maxDeg = max(degrees)
		strengths = colSums(monthAdj)
		st = sum(monthAdj)
		maxS = max(strengths)
		sparsity = (N^2-sum(monthAdj==0))/N^2

		## Entering stats
		monthIndx = 12*(i-1) + j
		monthMat[monthIndx,] = c(nZeroDegs,
						maxDeg,
						maxS,
						st,
						sparsity)
		
	}

	## Calc stats for yrAdj
	degrees = colSums(yrAdj>0)
	nZeroDegs = sum(degrees==0)
	maxDeg = max(degrees)
	strengths = colSums(yrAdj)
	st = sum(yrAdj)
	maxS = max(strengths)
	sparsity = (N^2-sum(yrAdj==0))/N^2

	## Entering stats for yrAdj
	yearMat[i,] = c(nZeroDegs,
				maxDeg,
				maxS,
				st,
				sparsity)

}

save(monthMat,yearMat,firstYear,lastYear,file=file.path(airDir,"dataInfo.RData"))
