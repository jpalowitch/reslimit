library(grDevices)

# Getting plot-able airports (from lat/long file)
code2LL = read.csv("http://stats.johnpalowitch.com/data-sets/CCME_analyses/airports/code2LL.csv",
                   header = TRUE, stringsAsFactors = FALSE)
prY <- range(code2LL$Latitude)
prX <- range(code2LL$Longitude)
pAblePorts <- code2LL$locationID

####################################
##	Makes dot plot of given codes
####################################
plotCodes = function(codes,cex = 0.5,pch = 16){
	mch = match(codes,code2LL$locationID)
	missing = which(is.na(mch))
	mch = mch[!is.na(mch)]
	points(code2LL$Longitude[mch],code2LL$Latitude[mch],pch=pch,cex=cex)
	cat(length(missing)," codes with no plot data\n")
}

####################################
##	Makes colored plot of partition
####################################

plotPart = function(partCodes,symbols = NULL,colPal = NULL,cex = 2){
	K = length(partCodes)
	sym_col_ord <- order(unlist(lapply(partCodes, length)), decreasing = TRUE)
	if(is.null(symbols))
		symbols = sym_col_ord
	if(is.null(colPal))
		colPal = rainbow(K)[sym_col_ord]


	for(i in 1:K){
		codes = partCodes[[i]]
		mch = match(codes,code2LL$locationID)
		missing = which(is.na(mch))
		mch = mch[!is.na(mch)]
		points(code2LL$Longitude[mch],code2LL$Latitude[mch],col = colPal[symbols[i]],pch=symbols[i],cex = cex)
	}
}


####################################
##	Calculating optimal png dims
####################################

pngSize = 800
stndUnit = prX[2]-prX[1]
asp = c(1,(prY[2]-prY[1])/stndUnit)
pngW = pngSize*asp[1]
pngH = pngSize*asp[2]


####################################
##	Wrapper function 
####################################

plotResults = function(output, codeList, foutTag, plotTitle, symbols = NULL, colPal = NULL, commCex = 2, bgCex = 0.5, bgPch = 16,filterSize = 1){

	commLengths = unlist(lapply(output$communities,length))

	### Filtering community sizes
	if(filterSize>1){
		tooSmall = output$communities[commLengths<filterSize]
		output$background = c(output$background,unlist(tooSmall))
		output$communities = output$communities[commLengths>=filterSize]
	}

	### Finding largest 25 communities
	if(length(output$communities)>25)
		output$communities = output$communities[order(commLengths,decreasing = TRUE)[1:25]]
	

	### Transforming numbers into codes
	output$communities = lapply(output$communities,function(indx)codeList[indx])

	png(paste0(foutTag,".png"),width = pngW, height = pngH)
	plot(prX[1]*1.5,prY[1]*1.5,col="white",
	main = plotTitle,
      xlim = rev(prX),
      ylim = prY,
	xlab = "Longitutde",
	ylab = "Latitude")
	plotCodes(codeList,cex = bgCex,pch = bgPch)
    if(length(output$communities)>1)
	    plotPart(output$communities,symbols,colPal,cex = commCex)
	dev.off()

}
	