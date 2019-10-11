### RegionSelection [Region, RegionFile, Extent] (selecting region and extent from
### shapefiles)
RegionSelection <- function(Region, RegionFile, Extent) {
  ## loading shapefiles
  Shapes <- readOGR(Dir.Mask, "ne_50m_admin_0_countries", verbose = FALSE)
  ## selecting region from shapefile run global analysis read user-defined extent
  ## (if applicable)
  if (Region == "Global") {
    if (is.null(Extent)) {
      area <- extent(-180, 180, -90, 90)
    } else {
      area <- Extent
    }
    location <- 1:length(Shapes) # selecting all countries contained within the shapefile
  } else {
    Where <- Region # countries to consider
    location <- NA # position vector in shapefile list
    for (i in 1:length(Where)) {
      # select region from Shapefiles
      location[i] <- which(as.vector(Shapes$NAME) == Where[i])
    }
    if (is.null(Extent)) {
      # read user-defined extent (if applicable)
      area <- extent(Shapes[location, ])
    } else {
      area <- Extent
    }
  }
  if (is.null(RegionFile)) {
    # if no file name has been specified
    RegionFile <- toString(Region) # take name of region
  }
  # returning parameters
  return(list(area, location, RegionFile))
}

### Fun_NamesRas [raster, ClimVar, ClimVar2]
# (assigning layer names to model rasters) ----
Fun_NamesRas <- function(raster, ClimVar, ClimVar2, rasiter = 1){
  names(raster) <- c(paste("Most informative", ClimVar[[rasiter]],"lag", sep=" "),
                     "Model AICs", "Model p-value", "Antecedent NDVI (c_NDVI)",
                     paste("Antecedent" , ClimVar[[rasiter]], "(c_clim)", sep=" "),
                     paste("Antecedent" , ClimVar2[[rasiter]], "(c_clim2)", sep=" "),
                     "Explained Variance", "Variance (NDVI)", "Variance (Shared)",
                     paste("Variance (" , ClimVar[[rasiter]], ")", sep=""))
  return(raster)} # Fun_NamesRas end