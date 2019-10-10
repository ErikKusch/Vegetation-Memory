setwd(Dir.Gimms)
gimms_files <- updateInventory
####--------------- RasterGIMMs [from, to] ----
RasterGIMMs <- function(from, to){
  print("#################################################")
  print(paste("Rasterising GIMMs NDVI data from ", from, " to ", to, sep=""))
  setwd(Dir.Gimms) # set working directory to base GIMMs folder
  # PREPARING DATA----
  invisible(capture.output(gimms_files <- downloadGimms(x = as.Date(paste(from,"-01-01",sep="")), # start date
                                                        y = as.Date(paste(to,"-12-31",sep="")), # end date
                                                        dsn = Dir.Gimms))) # where to store the files
  gimms_raster <- rasterizeGimms(x = gimms_files, remove_header = TRUE) # rasterising
  indices <- monthlyIndices(gimms_files) # extract month indices from file list (should be two of each)
  gimms_raster_mvc <- monthlyComposite(gimms_raster, indices = indices) # create composites according to indices
  # Fix NDVI misbehaviours
  gimms_raster_mvc[gimms_raster_mvc<0] <- 0 # set threshold for barren land (NDVI<0)
  gimms_raster_mvc[gimms_raster_mvc>1] <- 1 # set threshold for saturated NDVI (NDVI > 1)
  # Indices
  Years <- rep(seq(from, to, 1), each = 12) # The year corresponding to each month in the stack
  names(gimms_raster_mvc) <- paste(month.abb, Years, sep = "") # create names for rasters
  # SAVING DATA----
  writeRaster(gimms_raster_mvc, paste(Dir.Gimms.Monthly, "/GlobalNDVI_", from, to, sep=""),
              overwrite=TRUE, format="CDF", varname="GIMMsNDVI")
  setwd(mainDir)}# end of RasterGIMMS-function
####--------------- CombineCDFS [Region, RegionFile, Extent] ----
CombineCDFs <- function(Region, RegionFile, Extent){
  print("#################################################")
  print(paste("Producing cropped GIMMs NDVI raster stacks for ", RegionFile, sep=""))
  # SELECTING FILES----
  files <- list.files(Dir.Gimms.Monthly)
  files.pos <- grep("GlobalNDVI", files)
  # REGION SELECTION----
  Shapes <- readOGR(Dir.Mask,'ne_50m_admin_0_countries', verbose = FALSE)
  RegObj <- RegionSelection(Region = Region, RegionFile = RegionFile, Extent = Extent)
  area <- RegObj[[1]]
  location <- RegObj[[2]]
  RegionFile <- RegObj[[3]]
  # LOADING, CROPPING AND MASKING----
  setwd(Dir.Gimms.Monthly)
  ras <- list()
  for(i in 1:length(files.pos)){
    rasinter <- brick(files[files.pos[i]]) # load i-th ndvi file
    rasinter <- crop(rasinter, area) # cropping to extent
    rasinter <- mask(rasinter, Shapes[location,]) # masking via Shapefile
    ras[[i]] <- rasinter # save masked ndvi to list
  }
  ras <- brick(ras) # create one big regional ndvi raster
  # SAVING DATA----
  writeRaster(ras, paste(Dir.Gimms.Monthly, "/NDVI_", RegionFile, sep=""),
              overwrite=TRUE, format="CDF", varname="NDVI",
              longname= paste("Monthly NDVI means across ", Region, sep=""))
  setwd(mainDir)} # CombineCDFs