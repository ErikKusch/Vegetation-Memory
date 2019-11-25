# Fun_PFTs(Regions = list(c("Portugal", "Spain", "France", "Andorra"), "Brazil", "Australia"),
#          RegionFiles = list("Iberian Region", "Caatinga", "Australia"),
#          Extents = list(extent(-10,10,35,52), extent(-50,-34,-23,0), NULL),
#          From = 1982, To = 2015, Occ = "Download")


# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# print("CALCULATING SPECIES SPECIFIC-TRAIT MEANS")
# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# # calculating species-specific trait means
# if (!file.exists(paste(Dir.TRY, "/SpeciesTraits.RData", sep = ""))) {
#   PFTs() # species-specific trait means
# } else {
#   print("Species-specific trait means already calculated")
# }
# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# print("OBTAINING SPECIES OCCURENCE RECORDS FROM GBIF")
# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# # figuring out ISO3166 codes for selected regions
# ISO3166_df <- read.csv(paste(Dir.Mask, "/ISO3166.csv", sep = "")) # read ISO3166 country code list
# CountCodes <- ISO3166_df$ISO.3166.ALPHA.2[which(ISO3166_df$Country %in% unlist(Regions))]
# CountCodes <- toString(CountCodes)
# CountCodes <- gsub(pattern = ", ", replacement = ";", x = CountCodes)
# # Download occurence files from GBIF
# if (Occ == "Download") {
#   GbifStat <- NULL
#   # species occurences
#   Attempt <- 0
#   while (is.null(GbifStat) || GbifStat != "Done") {
#     if (Attempt > 0) {
#       print("Encountered an error and starting the downloading of GBIF occurence data again.
#               This is usually due to issues with the GBIF server connection and you don't have
#               to worry as long as your internet connection is stable.")
#     }
#     Attempt <- Attempt + 1
#     try(GbifStat <- DistMaps(Species = "All", Years = From:To, CountCodes = CountCodes))
#   }
# } else {
#   print("Occurence data will not be downloaded according to function call.")
# }
# # generate mean rasters of PFTs across specifiec regions looping over regions
# for (CompReg in 1:length(Regions)) {
#   if (!file.exists(paste(Dir.TRY, "/TRY-", RegionFiles[[CompReg]], ".nc", sep = ""))) {
#     PFTRasters(Region = Regions[[CompReg]], RegionFile = RegionFiles[[CompReg]],
#                Extent = Extents[[CompReg]], CountCodes = CountCodes)
#   } else {
#     print(paste("TRY data already aggregated to mean raster for region ",
#                 RegionFiles[[CompReg]], sep = ""))
#   }
# } # region-loop








####--------------- PFTs [] 
# (loading data, building species-specific trait means and saving the result) ----
PFTs <- function(){
  # LOADING DATA ----
  ## NDVI (reference raster) ---
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/GlobalNDVI_20112015.nc", sep=""))
  ref_ras <- NDVI_ras[[6]]
  ref_ras[which(values(ref_ras) > -1)] <- 8888 # identify land pixels
  ## Master PFT data from TRY
  PFT_Master <- read.table(file = paste(Dir.TRY, "/4704.txt", sep=""), stringsAsFactors = FALSE, fill = TRUE,
                           sep="\t", header = TRUE)
  ## Extracting necessary data to handle smaller data frame
  PFTs_df <- data.frame(Species = PFT_Master$SpeciesName, ObsID = PFT_Master$ObservationID,
                        Variable = PFT_Master$DataName, Value = PFT_Master$OrigValueStr, Unit = PFT_Master$UnitName)
  ## fixing factor to numeric
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  PFTs_df$Value <- as.numeric.factor(PFTs_df$Value)
  # CALCULATING SPECIES-SPECIFIC MEAN TRAIT VALUES ----
  ## Preparations for calculations
  Species <- unique(PFTs_df$Species) # all species to consider
  FullSpec_df <- data.frame(Species = NA, Nmass = NA, Height = NA)
  ## Calculations
  for(Iter in 1:length(Species)){ # loop over all species that need consideration
    T_Begin <- Sys.time() # read start time (needed for expected finishing time)
    ### Nitrogen mean
    Nitro <- which(PFTs_df$Species == Species[Iter] & # positions of nitrogen rows for species
                     PFTs_df$Variable == "Leaf nitrogen content per dry mass (Nmass)")
    Mean_Nitro <- mean(PFTs_df$Value[Nitro], na.rm = TRUE) # mean
    ### Height mean
    Height <- which(PFTs_df$Species == Species[Iter] & # positions of height rows for species
                      PFTs_df$Variable == "Plant height vegetative")
    Mean_Height <- mean(PFTs_df$Value[Height], na.rm = TRUE) # mean
    ### Combining data into a data frame
    Spec_df <- data.frame(Species = Species[Iter], Nmass = Mean_Nitro, Height = Mean_Height)
    FullSpec_df <- rbind(FullSpec_df, Spec_df)
    ### Updating progress bar
    if(Iter == 1){ # estimate on finishing time on first loop
      T_End <- Sys.time() # read end time
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # duration between date points
      ### Output to console
      print(paste("Calculating species-specific trait mean values should finish around ",
                  as.POSIXlt(T_Begin + Duration*length(Species), tz = Sys.timezone(location=TRUE)), sep=""))
      ### Removing empty initial row
      FullSpec_df <- FullSpec_df[-1,]}# end of estimator if-statement
    pb <- txtProgressBar(min = 0, max = length(Species), style = 3) # Setting up a progress bar
    setTxtProgressBar(pb, Iter)}# end of species-specific mean trait value calculation for-loop
  FullNA <- which(is.na(FullSpec_df$Nmass) & is.na(FullSpec_df$Height)) # identifying species where both means are NA
  CorrectedSpec_df <- FullSpec_df[-FullNA,] # removing full NA species records
  # Save the species-specific and NA-freed data
  save(CorrectedSpec_df, file = paste(Dir.TRY, "/SpeciesTraits.RData", sep=""))
  setwd(mainDir)
  # BUILDING RASTERS FROM RAW TRAIT MEASURES ----
  ## create empty data frame and filling it
  Locs_df <- data.frame(H = NA, Nmass = NA, Lat = NA, Lon = NA)
  `%nin%` = Negate(`%in%`) # create a 'not in' statement
  ## eliminate observations with less than three records (these can't have enough data)
  IDs <- PFTs_df$ObsID
  exclude <- names(table(IDs))[which(table(IDs) < 3)]
  IDs <- unique(IDs)[which(as.character(unique(IDs)) %nin% exclude)]
  ## loop over all observations with enough data
  for(i in 1:length(IDs)){
    T_Begin <- Sys.time() # read start time (needed for expected finishing time)
    Iter_df <- PFTs_df[which(PFTs_df$ObsID == IDs[i]),] # extract all data for current observation
    ## data extraction
    Nmass <- Iter_df$Value[which(Iter_df$Variable == "Leaf nitrogen content per dry mass (Nmass)")]
    if(length(Nmass) == 0){Nmass <- NA}
    H <- Iter_df$Value[which(Iter_df$Variable == "Plant height vegetative")]
    if(length(H) == 0){H <- NA}
    Lat <- as.numeric(Iter_df$Value[which(Iter_df$Variable == "Latitude")])
    Lon <- as.numeric(Iter_df$Value[which(Iter_df$Variable == "Longitude")])
    Locs_df <- rbind(Locs_df, c(H, Nmass, Lat, Lon)) # bin data
    ## Updating progress bar
    if(i == 1){ # estimate on finishing time on first loop
      T_End <- Sys.time() # read end time
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # duration between date points
      ### Output to console
      print(paste("Extracting raw geo-referenced data should finish around ",
                  as.POSIXlt(T_Begin + Duration*length(IDs), tz = Sys.timezone(location=TRUE)), sep=""))
      ### Removing empty initial row
      Locs_df <- Locs_df[-1,]
      pb <- txtProgressBar(min = 0, max = length(IDs), style = 3)}
    setTxtProgressBar(pb, i)
  } # obvservation loop
  ## converting to spatial points objects
  H_df <- na.omit(Locs_df[which(!is.na(Locs_df$H)),-2])
  H_pts <- data.frame(y = H_df$Lat, x = H_df$Lon, z = H_df$H)
  H_pts <- na.omit(H_pts)
  coordinates(H_pts) = ~x+y # convert x and y to coordinates
  Nmass_df <- na.omit(Locs_df[which(!is.na(Locs_df$Nmass)),-1])
  Nmass_pts <- data.frame(y = Nmass_df$Lat, x = Nmass_df$Lon, z = Nmass_df$Nmass)
  Nmass_pts <- na.omit(Nmass_pts)
  coordinates(Nmass_pts) = ~x+y # convert x and y to coordinates
  ## rasterising
  rast <- raster(ext=extent(ref_ras), resolution=res(ref_ras)) # create raster to be filled
  H_rasOut <- rasterize(x = H_pts, y = rast, field = H_pts$z, fun = mean) # rasterize irregular points
  Nmass_rasOut <- rasterize(x = Nmass_pts, y = rast, field = Nmass_pts$z, fun = mean) # rasterize irregular points
  Means_ras <- brick(H_rasOut, Nmass_rasOut)
  ## saving data
  writeRaster(x=Means_ras, filename = paste(Dir.TRY,"/RawTRY-Global", sep=""),overwrite=TRUE, format="CDF")
}# end of PFTs-function

####--------------- DistMaps [Species, Extent, Years, CountCodes]
# (Obtaining occurence data via GBIF, rasterising, saving the raster, limitting to a region) ----
DistMaps <- function(Species, Years, CountCodes){
  print("#################################################")
  print(paste("Downloading occurence data of species: ", Species, " across: ", CountCodes, sep=""))
  # LOADING SPECIES DATA FRAME ----
  load(paste(Dir.TRY, "/SpeciesTraits.RData", sep="")) # load data frame 'CorrectedSpec_df'
  if(Species == "All"){ # selecting all species contained in 'CorrectedSpec_df'
    Species <- sort(CorrectedSpec_df$Species)
    Species <- Species[-which(Species == "-")] # remove this error of a species name
    SP <- "All"
  }else{
    SP <- "Dummy" # only used if single species is targeted
  }
  # GLOBAL REFERENCE DATA (needed for rasterising and masking) ----
  ## NDVI (reference raster) ---
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/GlobalNDVI_20112015.nc", sep=""))
  ref_ras <- NDVI_ras[[6]]
  ref_ras[which(values(ref_ras) > -1)] <- 8888 # identify land pixels
  # IDENTIFYING GBIF KEY(S) ----
  # if data is already present and all species are sought-after
  if("SpeciesGBIFKeys.rda" %in% list.files(path=Dir.TRY) & SP == "All"){
    print("Loading species-specific GBIF keys from local storage") # output to console
    load(paste(Dir.TRY, "/SpeciesGBIFKeys.rda", sep=""))
  }else{
    ## Preparations ---
    Key_vec <- NA # create empty vector for gbif key(s)
    Species_Pres <- NA # create empty vector for all species which we have occurence records for
    print("Identifying species-specific GBIF keys from GBIF repository") # output to console
    pb <- txtProgressBar(min = 0, max = length(Species), style = 3) # Setting up a progress bar
    for(Iter in 1:length(Species)){ # cycle through all species specified
      key <- name_suggest(q=Species[Iter], rank='species')$key[1] # pull gbif key
      Key_vec <- c(Key_vec, key) # append key to key vector
      if(!is.null(key)){ # if we have occurence data in the gbif records
        Species_Pres <- c(Species_Pres, Species[Iter]) # append species to species vector
      }
      setTxtProgressBar(pb, Iter)}
    ## Fixing vectors
    Key_vec <- Key_vec[-1] # Removing empty initial element
    Species_Pres <- Species_Pres[-1] # removing empty initial element
    # Save the species-specific and NA-freed data
    if(SP == "All"){ # if we want all species, we might as well save the names and key objects for later saving of time
      save(list = c("Species_Pres", "Key_vec"), file = paste(Dir.TRY, "/SpeciesGBIFKeys.rda", sep=""))}} # GBIF keys
  # OCCURENCE DATA ----
  ## Preparation
  print("Downloading species-specific occurence records from GBIF") # output to console
  pb <- txtProgressBar(min = 0, max = length(Species_Pres), style = 3) # Setting up a progress bar
  ## If an error occured previously
  if("Breakage.txt" %in% list.files(path = mainDir)){ # this file is only present if the run finished prematurely
    OccIter <- read.table(paste(mainDir, "/Breakage.txt", sep=""))[1,1] # position at which it failed previously - 1
  }else{
    OccIter <- 0 # set to 0 if it didn't fail previously
  }
  if(OccIter > 1){ # if previous run (OccIter) failed at the second step or later
    Start <- OccIter + 1 # start from where it failed, OccIter is the last one that got done
  }else{ # if it failed at the first one
    Start <- 1 # start at the first species
  }
  for(OccIter in Start:length(Key_vec)){ # cycling through all species to obtain occurence data
    # if species name cannot be put into a file name due to special characters, this excludes 61 species records
    if(grepl('[^[:alnum:]. .).(.-]', Species_Pres[OccIter]) == TRUE){
      next()}
    if(paste(Species_Pres[OccIter], "_", CountCodes, ".rda", sep="") %in% list.files(path=Dir.OCCs)){
    }else{ # data not present locally yet
      ## Downloading Data
      key <- Key_vec[OccIter] # select GBIF key
      Gbif <- occ_data(key, limit=200000, hasCoordinate = TRUE, year = Years,
                       hasGeospatialIssue = FALSE, country = CountCodes) # download data
      ## Dealing with separate data frames of years
      BaseOcc <- rep(NA, 3) # create empty vector
      BaseOcc_df <- t(as.data.frame(BaseOcc)) # make empty vector into empty data frame
      colnames(BaseOcc_df) <- c("decimalLatitude", "decimalLongitude", "year") # set column names
      for(i in 1:length(Years)){
        # create a data frame of latitude and longitude records of currently iterated year
        GbifFrame <- data.frame(decimalLatitude = Gbif[[i]]$data$decimalLatitude,
                                decimalLongitude = Gbif[[i]]$data$decimalLongitude,
                                year = rep(Years[i], length(Gbif[[i]]$data$decimalLatitude)))
        BaseOcc_df <- rbind(BaseOcc_df, GbifFrame)}
      ## Sanity check
      if(dim(BaseOcc_df)[1] == 1){ # if there is no occurence data
        next()}
      BaseOcc_df <- BaseOcc_df[-1,] # remove initial NA row
      ## Saving data frame
      save(BaseOcc_df, file = paste(Dir.OCCs,"/", Species_Pres[[OccIter]], "_", CountCodes,".rda",sep=""))}
    setTxtProgressBar(pb, OccIter) # pdate progress bar
    # save current iteration number to disk (used for jumping right back in if errors occur)
    write.table(OccIter, file = paste(mainDir, "/Breakage.txt", sep=""))
  }# occurence data loop
  file.remove(paste(mainDir, "/Breakage.txt", sep=""))
  setwd(mainDir)
  GbifStat <- "Done"
  return(GbifStat)
}# end of Mapping function

####--------------- PFTRasters [Region, Extent, RegionFile, CountCodes]
# (loading data, building species-specific trait mean rasters for study regions) ----
PFTRasters <- function(Region, Extent, RegionFile, CountCodes){
  print("#################################################")
  print(paste("Building mean trait rasters across ", RegionFile,sep=""))
  load(paste(Dir.TRY, "/SpeciesTraits.RData", sep="")) # load data
  RawTry_ras <- brick(paste(Dir.TRY,"/RawTRY-Global.nc", sep=""))
  # GLOBAL REFERENCE DATA (needed for rasterising and masking) ----
  ## NDVI (reference raster) ---
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/GlobalNDVI_20112015.nc", sep=""))
  ref_ras <- NDVI_ras[[6]]
  ref_ras[which(values(ref_ras) > -1)] <- 8888 # identify land pixels
  # REGION SELECTION----
  Shapes <- readOGR(Dir.Mask,'ne_50m_admin_0_countries', verbose = FALSE)
  RegObj <- RegionSelection(Region = Region, RegionFile = RegionFile, Extent = Extent)
  area <- RegObj[[1]]
  location <- RegObj[[2]]
  RegionFile <- RegObj[[3]]
  # CROPPING AND MASKING ----
  ## Reference cropping and masking
  ref_rasC <- crop(ref_ras, area) # cropping to extent
  ref_rasF <- mask(ref_rasC, Shapes[location,]) # masking via Shapefile
  # RAW TRY DATA ----
  RawTry_rasC <- crop(RawTry_ras, area) # cropping to extent
  RawTry_rasF <- mask(RawTry_rasC, Shapes[location,]) # masking via Shapefile
  writeRaster(x=RawTry_rasF, filename = paste(Dir.TRY,"/RawTRY-",RegionFile, sep=""),overwrite=TRUE, format="CDF")
  # CALCULATING MEAN RASTERS WITH DISTRIBUTION MAPS ----
  # create empty mean raster
  BaseMeans <- ref_rasF
  values(BaseMeans)[!is.na(values(BaseMeans))] <- 0
  # build brick for mean calculations
  BaseMeans <- brick(BaseMeans, BaseMeans, BaseMeans, BaseMeans)
  names(BaseMeans) <- c("Height", "NMass", "HCount", "NCount")
  # progress bar
  pb <- txtProgressBar(min = 0, max = length(list.files(Dir.OCCs)), style = 3) #
  # looping over all .rda occurence files previously downloaded
  for(OccRast in 1:length(list.files(Dir.OCCs))){
    # OCCURENCE ----
    load(paste(Dir.OCCs, "/", list.files(Dir.OCCs)[OccRast], sep=""))
    ## Converting to SpatialPoints
    pts <- data.frame(y = BaseOcc_df$decimalLatitude, x = BaseOcc_df$decimalLongitude,
                      z = rep(1, length(BaseOcc_df$decimalLongitude)))
    pts <- na.omit(pts) # remove NA rows
    coordinates(pts) = ~x+y # convert x and y to coordinates
    # RASTERISING ----
    # create raster to be filled
    rast <- raster(ext=extent(ref_ras), resolution=res(ref_ras))
    # rasterize irregular points
    # we use a mean function here to regularly grid the irregular input points
    rasOut<-rasterize(x = pts, y = rast, field = pts$z, fun = max)
    ## Occurence cropping and masking
    rasC <- crop(rasOut, area) # cropping to extent
    rasF <- mask(rasC, Shapes[location,]) # masking via Shapefile
    # TRAIT MEANS ----
    # loading data of currently iterated on species
    Grep <- list.files(Dir.OCCs)[OccRast]
    Grep <- gsub(x = Grep, pattern = CountCodes, replacement = "")
    Grep <- gsub(x = Grep, pattern = "_.rda", replacement = "")
    
    NMass <- CorrectedSpec_df$Nmass[which(CorrectedSpec_df$Species == Grep)]
    Height <- CorrectedSpec_df$Height[which(CorrectedSpec_df$Species == Grep)]
    # BUILDING MAP
    Identifier <- which(!is.na(values(rasF)))
    if(length(NMass) != 0){
      if(!is.nan(NMass)){ # add current species-NMass to raster layer and bump up count by 1
        values(BaseMeans$NMass)[Identifier] <- values(BaseMeans$NMass)[Identifier] + NMass
        values(BaseMeans$NCount)[Identifier] <- values(BaseMeans$NCount)[Identifier] + 1}}
    if(length(Height) != 0){
      if(!is.nan(Height)){ # add current species-Height to raster layer and bump up count by 1
        values(BaseMeans$Height)[Identifier] <- values(BaseMeans$Height)[Identifier] + Height
        values(BaseMeans$HCount)[Identifier] <- values(BaseMeans$HCount)[Identifier] + 1}}
    setTxtProgressBar(pb, OccRast) # update progress bar
  } # OccRast-loop
  # CALCULATE MEANS ----
  TestHeight <- BaseMeans$Height/BaseMeans$HCount
  values(TestHeight)[which(values(TestHeight) > quantile(values(TestHeight), .95, na.rm = TRUE))] <- NA
  TestNMass <- BaseMeans$NMass/BaseMeans$NCount
  values(TestNMass)[which(values(TestNMass) > quantile(values(TestNMass), .95, na.rm = TRUE))] <- NA
  Means_ras <- brick(TestHeight, TestNMass)
  # SAVING DATA ----
  writeRaster(x=Means_ras, filename = paste(Dir.TRY,"/TRY-",RegionFile, sep=""),overwrite=TRUE, format="CDF")
}# PFTRasters