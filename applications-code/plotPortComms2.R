library(grDevices)
library(ggplot2)
library(igraph)
library(ggmap)
library(grid)

## Getting plot-able airports (from lat/long file)
code2LL = read.csv("http://stats.johnpalowitch.com/data-sets/CCME_analyses/airports/code2LL.csv",
                   header = TRUE, stringsAsFactors = FALSE)
lat_range <- range(code2LL$Latitude)
lon_range <- -range(code2LL$Longitude)
plot_center <- c(-mean(lon_range), mean(lat_range))
prY = range(code2LL$Latitude)
prX = range(code2LL$Longitude)
pAblePorts = code2LL$locationID
full_loc <- c(lon_range[2], lat_range[1], lon_range[1], lat_range[2])

# Define custom symbol order
sym_ord <- c(11, 20, 17, 7, 8, 9)
sym_ord <- c(sym_ord, setdiff(1:25, sym_ord))

# Distinct colors:
dcolors <- c(
"#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
  
  # Color choices:
  # plot(rep(1:13, each = 8), rep(1:8, 13), col = dcolors, pch = 15, cex = 3)
  # text(rep(1:13, each = 8), rep(1:8, 13), labels = as.character(1:104))
  cchoices <- c(2:5, 42, 78, 24, 6:7, 9:12, 14:15, 21, 28:31)
  
colPal <- dcolors[cchoices]

####################################
##	Wrapper function 
####################################

plotResults = function(output, codeList, strengths = NULL, adj_names = NULL,
                       plotTitle, filterSize = 1){

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

	
	p <- qmap(full_loc, maptype = "toner-lite", color = "bw") + 
	 geom_point(data = mship_df, 
	            aes(x = lon, 
	                y = lat, 
	                #size = str, 
	                size = 10,
	                color = col,
	                shape = sym)) +
	 scale_color_identity() + 
	 scale_shape_identity() + 
	 #scale_size(range=range(mship_df$str) / 2 , guide="none") +
	  scale_size_identity(guide="none") + 
	  xlab("Longitude") + 
	  ylab("Latitude") +
	  theme(axis.title.y = element_text(size = rel(3)),
	        axis.title.x = element_text(size = rel(3)),
	        axis.text.x = element_text(size = rel(3)),
	        axis.text.y = element_text(size = rel(3)))
	return(p)
	

}
	