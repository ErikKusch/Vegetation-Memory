####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions

####--------------- DATA ----
ERA <- raster(paste(Dir.ERA.Monthly, "Tair_mean_GlobalDrylands_011981_122015.nc", sep="/"))
NDVI <- raster(paste(Dir.Gimms.Monthly, "GlobalDrylandsNDVI.nc", sep="/"))
Back <- brick(paste(Dir.Gimms.Monthly, "GlobalNDVI_19821985.nc", sep="/"))[[6]]

####--------------- PLOTS ----
plot(Back, col="grey", colNA = "black", legend=FALSE, axes=TRUE, 
     main = "Air Temperature [K]")
plot(ERA, col=heat.colors(100), legend=TRUE, axes=FALSE, add=TRUE, 
     legend.width=2, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1)

plot(Back, col="grey", colNA = "black", legend=FALSE, axes=TRUE, 
     main = "Normalized Difference Vegetation Index (NDVI)")
plot(NDVI, col=colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000), 
     legend=TRUE, axes=FALSE, add=TRUE, 
     legend.width=2, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1)