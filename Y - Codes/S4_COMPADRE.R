####--------------- RasterCOMPADRE [Variable, Region, RegionFile, Extent]
# (Selecting specified COMPADRE data, rasterising, saving the raster, limitting to a region) ----
RasterCOMPADRE <- function(Variable, Region, RegionFile, Extent){
  print("#################################################")
  print(paste("Rasterising COMPADRE ", Variable, " across ", RegionFile, sep=""))
  # LOADING DATA ----
  Compadre_df <- read.csv(paste(Dir.Compadre, "/derived_demographicsAllCOMPADRE.csv", sep="")) # load data frame
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/GlobalNDVI_20112015.nc", sep="")) # reference raster
  ref_ras <- NDVI_ras[[6]] # select only one years data
  ref_ras[which(values(ref_ras) > -1)] <- 8888 # select only land pixels and set them to -8888
  # REGION SELECTION----
  Shapes <- readOGR(Dir.Mask,'ne_50m_admin_0_countries', verbose = FALSE)
  RegObj <- RegionSelection(Region = Region, RegionFile = RegionFile, Extent = Extent)
  area <- RegObj[[1]]
  location <- RegObj[[2]]
  RegionFile <- RegObj[[3]]
  # DATA MANIPULATION ----
  if(Variable == "FastSlow"){ # analysis of PCA axes
    FSVars <- c("GenT", "H", "La", "GrowSSD", "ShriSSD", "RepSSD", "S", "R0", "Lmean")
    VariableCol <- match(FSVars, colnames(Compadre_df))
    FSLoads <- list(c(.87,.53,.7,-.8,.04,-.81,-.25,-.03,.12), # PCA 1 according to Salguero-Gomez, 2017
                    c(.15,.27,.28,-.05,-.79,.32,.65,.7,.26)) # PCA 2 according to Salguero-Gomez, 2017
      pts <- na.omit(data.frame(y = Compadre_df$Lat, x = Compadre_df$Lon, z = Compadre_df[,VariableCol]))
    ## PCA calculations
    PCA1_df <- t(t(pts[,-1:-2]) * FSLoads[[1]]) # multiplying by first axis loadings
    PCA1_df <- rowSums(PCA1_df) # build sums for single index along PCA 1
    PCA1_df <- cbind(pts[,1:2], PCA1_df) # binding with coordinates
    coordinates(PCA1_df) = ~x+y # convert x and y to coordinates
    PCA2_df <- t(t(pts[,-1:-2]) * FSLoads[[2]]) # multiplying by second axis loadings
    PCA2_df <- rowSums(PCA2_df) # build sums for single index along PCA 2
    PCA2_df <- cbind(pts[,1:2], PCA2_df) # binding with coordinates
    coordinates(PCA2_df) = ~x+y # convert x and y to coordinates
    ## Rasterizing ----
    rast <- raster(ext=extent(ref_ras), resolution=res(ref_ras)) # create raster to be filled
    rasOut1 <-rasterize(x = PCA1_df, y = rast, field = PCA1_df$PCA1_df, fun = mean) # rasterize irregular points
    rasOut2 <-rasterize(x = PCA2_df, y = rast, field = PCA2_df$PCA2_df, fun = mean) # rasterize irregular points
    rasOut <- brick(rasOut1, rasOut2)
    names(rasOut) <- c("FS PCA1", "FS PCA2")
  }else{ # single variable desired
    VariableCol <- which(colnames(Compadre_df) == Variable) # figure out the position of the desired Variable
    pts <- na.omit(data.frame(y = Compadre_df$Lat, x = Compadre_df$Lon, z = Compadre_df[,VariableCol]))
    coordinates(pts) = ~x+y # convert x and y to coordinates
    rast <- raster(ext=extent(ref_ras), resolution=res(ref_ras)) # create raster to be filled
    rasOut<-rasterize(x = pts, y = rast, field = pts$z, fun = mean) # rasterize irregular points
  }
  # CROPPING AND MASKING ----
  rasC <- crop(rasOut, area) ## Occurence cropping and masking
  rasF <- mask(rasC, Shapes[location,]) # masking via Shapefile
  # DATA EXPORT ----
  Dir.Temp.Compadre <- paste(Dir.Compadre, Variable, sep="/")
  dir.create(Dir.Temp.Compadre)
  values(rasF)[which(values(rasF) == Inf)] <- NA # get rid off Inf values (when dealing with Rho)
  invisible(writeRaster(rasF, filename = paste(Dir.Temp.Compadre, "/",Variable,"_",RegionFile,sep=""),
                        overwrite=TRUE, format="CDF"))
}# end of RasterCOMPADRE

