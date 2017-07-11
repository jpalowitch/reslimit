########################################################################
## Function that:
##	- Takes a location of an OSLOM2 tp output
##	- Reads in the output 
##	- Returns the output in CCME format
########################################################################

osRead <- function (fin) {

	conn <- file(fin, open = "r")
	linn <- readLines(conn)
	osOut <- list(NULL)
	for (i in 1:length(linn)) {
	    readIn <- linn[i]
	    readIn <- strsplit(readIn," ")[[1]]
	    if (!((readIn[1] == "#module"))) {
	        readIn <- sapply(readIn, as.numeric)
	        names(readIn) <- NULL
	        osOut <- c(osOut, list(readIn))
	    }
	}
	close(conn)
	osOut <- osOut
	bgLocs <- unlist(lapply(osOut, length)) <= 1
	bg <- unlist(osOut[bgLocs])
	osOutFormat <- list("communities" = osOut[!bgLocs], "background" = bg)
	
	return(osOutFormat)

}