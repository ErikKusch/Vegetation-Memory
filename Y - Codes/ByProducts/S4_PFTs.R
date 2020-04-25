####--------------- BIEN [Trait, Region, RegionFile, Extent]
BIEN <- function(Trait, Region, RegionFile, Extent){
  print(paste("Handling", Trait, "data across", RegionFile))
  # PRELIMINARY ----
  # United States of America are treated as United States in BIEN data base, hence fixing needed
  if("United States of America" %in% Region){
    BIENReg <- Region
    BIENReg[which(Region == "United States of America")] <- "United States"
  }else{
    BIENReg <- Region
  }
  # LOADING DATA ----
  ## NDVI (reference raster) ---
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/GlobalNDVI_20112015.nc", sep=""))
  ref_ras <- NDVI_ras[[6]]
  ref_ras[which(values(ref_ras) > -1)] <- 8888 # identify land pixels
  ## BIEN data ---
  if(file.exists(paste(Dir.TRY, "/BIEN_", Trait, ".rds", sep=""))){
    BIEN_df <- readRDS(paste(Dir.TRY, "/BIEN_", Trait, ".rds", sep=""))
  }else{
    BIEN_df <- BIEN_trait_trait(trait = Trait, political.boundaries = TRUE, all.taxonomy = TRUE)
    saveRDS(BIEN_df, file = paste(Dir.TRY, "/BIEN_", Trait, ".rds", sep=""))
  }
  # DATA MANIPULATION ----
  Spec_BIEN <- BIEN_df[which(BIEN_df$country %in% BIENReg),] # limiting to Region (by country name)
  ## transforming into raster
  BIENpts <- na.omit(data.frame(y = Spec_BIEN$latitude, x = Spec_BIEN$longitude, 
                                z =as.numeric(Spec_BIEN$trait_value)))
  coordinates(BIENpts) = ~x+y # convert x and y to coordinates
  rast <- raster(ext=extent(ref_ras), resolution=res(ref_ras)) # create raster to be filled
  rasOut <- rasterize(x = BIENpts, y = rast, field = BIENpts$z, fun = mean) # rasterize irregular points
  # REGION SELECTION ----
  Shapes <- readOGR(Dir.Mask,'ne_50m_admin_0_countries', verbose = FALSE)
  RegObj <- RegionSelection(Region = Region, RegionFile = RegionFile, Extent = Extent)
  area <- RegObj[[1]]
  location <- RegObj[[2]]
  RegionFile <- RegObj[[3]]
  # CROPPING AND MASKING ----
  rasC <- crop(rasOut, area) ## Occurence cropping and masking
  rasF <- mask(rasC, Shapes[location,]) # masking via Shapefile
  # DATA SAVING ----
  invisible(writeRaster(rasF, filename = paste(Dir.TRY, "/", RegionFile, "_", Trait,sep=""),
                        overwrite=TRUE, format="CDF")) 
} # end of BIEN function