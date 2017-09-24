load(file.path("applications-results", "airports", "data", "adjList_cleaned.RData"))
load(file.path("applications-results", "airports", "data", "dataInfo.RData"))
foutTag <- file.path(curr_dir, yearNum, "test")
output <- list(communities = res$results,
               background = res$background)
filterSize <- 1
cex.pnt <- 3
methname <- "Recursive Modularity"

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
  labs(title = paste0(methname, ", ", toupper(fn), " ", yearNum)) + 
  theme(title = element_text(size = 40))

# Printing map and saving
png(paste0(foutTag,".png"), height = 1200, width = 1200)
print(pmap)
dev.off()