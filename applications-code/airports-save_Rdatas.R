# Reads in adjMats and
# makes .Rdatas for each
# one in correct folders

saveDir <- "applications-results/airports"

load("applications-results/airports/data/adjList_cleaned.RData")
load("applications-results/airports/data/dataInfo.RData")


## Extracting years
startYear <- as.numeric(names(adjList)[1])
endYear <- as.numeric(names(adjList)[length(adjList)])
yL <- endYear - startYear + 1

## Making month char list
months <- c("jan", "feb", "mar", "apr",
            "may", "jun", "jul", "aug",
            "sep", "oct", "nov", "dec")

for (y in 1:yL) {
    
    yearNum <- startYear + y - 1
    nMonths <- length(adjList[[y]]$months)
    
    for (m in 1:nMonths) {
        data <- adjList[[y]]$months[[m]]
        ysaveDir <- file.path(saveDir, yearNum)
        if (!dir.exists(ysaveDir))
          dir.create(ysaveDir, recursive = TRUE)
        fn <- file.path(ysaveDir, paste0(months[m], ".RData"))
        save(data, file = fn)
    }
    
    data <- adjList[[y]]$year
    fn <- file.path(ysaveDir, paste0("year.RData"))
    save(data, file = fn)
    
    
}

