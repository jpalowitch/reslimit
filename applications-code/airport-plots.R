source("applications-code/plotPortComms2.R")
rootDir = "applications-results/airports"
dataDir = file.path(rootDir, "data")
saveDir = file.path(rootDir, "plots")

if (!dir.exists(saveDir))
  dir.create(saveDir, recursive = TRUE)
monthNames = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
adj_options <- c(monthNames, "year")
fnList = c(monthNames,"year")
load(file.path(dataDir, "adjList_cleaned.RData"))

### Specify if you want bw plots or not
bw <- FALSE

### Specify whether or not to include plot titles
plotTitles = FALSE

### Set filter size (i.e. minimum cluster size)
filtSize = 1

### Set plot point size
cex.pnt <- 7

## Specify years
load(file.path(dataDir, "dataInfo.RData"))
startYear = firstYear
endYear = lastYear
yL = endYear - startYear + 1

  
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for (y in 1:yL) {
  
  cat("y <- ", y, "\n")
  
  yearNum <- startYear + y - 1
  allFiles <- list.files(file.path(rootDir, yearNum))
  outputFiles <- allFiles[grepl("output", allFiles)]
  outputFiles <- outputFiles[grepl(".RData", outputFiles)]
  fnList_y <- c(monthNames[1:length(adjList[[y]]$months)], "year")

  if (!dir.exists(file.path(saveDir, yearNum)))
    dir.create(file.path(saveDir, yearNum), recursive = TRUE)
  
	
	for (fn in fnList_y) {
	  
	  load(file.path(rootDir, yearNum, paste0(fn, ".RData")))
	  
	  # Getting fn-specific output files
	  fn_outputFiles <- outputFiles[grepl(paste0("_", fn, "."), 
	                                      outputFiles, fixed = TRUE)]
	  
	  # Getting codeList
	  codeList <- rownames(data)
	  
	  # Getting strengths and adjnames
	  if (fn == "year") {
	    strengths0 <- colSums(adjList[[y]]$year)
	    adj_names <- rownames(adjList[[y]]$year)
	  } else {
	    strengths0 <- colSums(adjList[[y]]$months[[which(monthNames == fn)]])
	    adj_names <- rownames(adjList[[y]]$months[[which(monthNames == fn)]])
	  }

		for (i in 1:length(fn_outputFiles)) {
            
		  
		  cat('i <- ', i, '\n')
			fname <- fn_outputFiles[i]
			load(file.path(rootDir, yearNum, fname))
      methname <- strsplit(fname, "_")[[1]][2]
      methname <- toupper(methname)
      if (methname == "SLPA")
        methname <- paste0(methname, "w")
            
      if(plotTitles){
			  plotTitle <- paste(fn, toString(yearNum))
      } else {
        plotTitle = ""
      }
			
			foutTag <- file.path(saveDir, yearNum, substr(fname, 1, nchar(fname) - 6))
			strengths  <- log10(strengths0) + 1
			output <- results
			filterSize <- 1
			
			commLengths = unlist(lapply(output$communities,length))
			
			### Filtering community sizes
			if(filterSize>1){
			  tooSmall = output$communities[commLengths<filterSize]
			  output$background = c(output$background,unlist(tooSmall))
			  output$communities = output$communities[commLengths>=filterSize]
			}
			
			### Finding largest 20 communities
			output$communities = output$communities[order(commLengths,decreasing = TRUE)[1:20]]
			
			
			### Transforming numbers into codes
			output$communities = lapply(output$communities,function(indx)codeList[indx])
			
			# Creating membership data frame
			mships <- lapply(1:length(output$communities),
			                 function (indx) rep(indx, length(output$communities[[indx]])))
			
			mship_df <- data.frame("code" = unlist(output$communities),
			                       "comm" = unlist(mships),
			                       stringsAsFactors = FALSE)
			
			# Removing those without L/L info
			LL_match <- match(mship_df$code, code2LL$locationID)
			mship_df <- subset(mship_df, !is.na(LL_match))
			
			# Adding flight info
			if (!is.null(strengths)) {
			  mship_df <- data.frame(mship_df,
			                         "str" = strengths[match(mship_df$code, adj_names)])
			} else {
			  mship_df <- data.frame(mship_df,
			                         "str" = rep(1, nrow(mship_df)))
			}
			
			# Adding LL info
			LL_match <- match(mship_df$code, code2LL$locationID)
			mship_df <- data.frame(mship_df,
			                       "lon" = -code2LL[LL_match, 3],
			                       "lat" = code2LL[LL_match, 2])
			
			
			# Inducing symbol order
			mship_df <- data.frame(mship_df,
			                       "sym" = as.integer(sym_ord[mship_df$comm]))
			
			# Inducing color choice
			#mship_df <- data.frame(mship_df, "col" = colPal[mship_df$comm])
			ggcolpal <- gg_color_hue(length(unique(mship_df$comm)))
			mship_df <- data.frame(mship_df, "col" = ggcolpal[mship_df$comm])
			
			# Randomize plotting
			mship_df <- mship_df[sample(nrow(mship_df)), ]
			
			# Factorizing community
			mship_df$comm <- factor(mship_df$comm)	
			
			# Making map
			mymap <- get_map(full_loc, maptype = "toner-lite")
			pmap <- ggmap(mymap, extent = "device")
			if (bw) {
			  pmap <- pmap + geom_point(aes(x = lon, y = lat, colour = "black", shape = sym),
			                            data = mship_df, size = cex.pnt) 
			} else {
			  pmap <- pmap + geom_point(aes(x = lon, y = lat, colour = col, shape = sym),
			                            data = mship_df, size = cex.pnt) 
			 }
			pmap <- pmap + guides(size = FALSE, colour = FALSE) + 
			  labs(title = "my title") + scale_colour_identity() + scale_shape_identity() +
			  labs(title = paste0(methname, ", ", toupper(fn), " ", seq(startYear, endYear)[y])) + 
			  theme(title = element_text(size = 40))
			
      # Printing map and saving
			png(paste0(foutTag,".png"), height = 1200, width = 1200)
			print(pmap)
			dev.off()
            
      rm(results, plotTitle)
 
		}
			
	}

}

sink()
