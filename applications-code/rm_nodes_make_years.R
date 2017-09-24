# Removes nodes with no connections --------------------------------------------
# Makes yearly matrices
# Saves output in /data

load("applications-results/airports/data/adjData.RData")
load("applications-results/airports/data/dataInfo.RData")

adjListNew <- adjList

# Extracting years
startYear <- as.numeric(names(adjList)[1])
endYear <- as.numeric(names(adjList)[length(adjList)])
yL <- endYear - startYear + 1

# Making month char list
for (y in 1:yL) {

	yearNum <- startYear + y - 1
	nMonths <- length(adjList[[y]]$months)
	yearMat <- matrix(0, nrow(adjList[[y]]$months[[1]]),
	                  ncol(adjList[[y]]$months[[1]]))
	rownames(yearMat) <- rownames(adjList[[y]]$months[[1]])
	colnames(yearMat) <- colnames(adjList[[y]]$months[[1]])
	
	for (m in 1:nMonths) {

		yearMat <- yearMat + adjList[[y]]$months[[m]]
	
		adjMat <- adjList[[y]]$months[[m]]
		posz <- colSums(adjMat) > 0
		adjMat <- adjMat[posz, posz]
		adjListNew[[y]]$months[[m]] <- adjMat
	}

	posz_y <- colSums(yearMat) > 0
	yearMat <- yearMat[posz_y, posz_y]
	adjListNew[[y]]$year <- yearMat

}

adjList <- adjListNew

save(adjList, file = "applications-results/airports/data/adjList_cleaned.RData")
	