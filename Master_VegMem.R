rm(list = ls()) # clearing environment
####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions
####--------------- VARIABLE VECTORS ----
ModVars <- c("Tair", "Qsoil1", "Qsoil2", "Qsoil3", "Qsoil4")
ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
###---------------- FUNCTIONS ---------------------------------------------------------
####--------------- Fun_Vegetation [Regions, RegionFiles, Extents, From, To, Lags, Cores]
# (selecting and preparing data, and calculating vegetation memory) ----
Fun_Vegetation <- function(Regions, RegionFiles, Extents, From, To, Lags, Cores) {
  FromY <- (From - ceiling(1/12 * max(Lags))) # Figuring out real start year after factoring in lags
  ##### GIMMs NDVI -----
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("HANDLING GIMMs NDVI DATA")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  setwd(mainDir)
  source(paste(Dir.Codes, "S1_GIMMs.R", sep="/"))
  # Download NDVI data for full years (from, to), turn it into monthly composite
  # rasters checking if all files are already there and not running this function
  # if they are
  froms <- c(1982, 1986, 1991, 1996, 2001, 2006, 2011)
  tos <- c(1985, 1990, 1995, 2000, 2005, 2010, 2015)
  for(RasGimms in 1:length(tos)){
    if(file.exists(paste(Dir.Gimms.Monthly, "/GlobalNDVI_", froms[RasGimms], tos[RasGimms], ".nc", sep=""))){
      print(paste("Global NDVI raster from", froms[RasGimms], "to", tos[RasGimms], "has already been established."))
      next()
    }else{
      RasterGIMMs(from = froms[RasGimms], to = tos[RasGimms])
    }
  }
  # Load composite NDVI data and limit to extent of a study region saving the
  # resulting data
  for (CombineRun in 1:length(Regions)) {
    # Checking if this particular data has been produced already
    if (file.exists(paste(Dir.Gimms.Monthly, "/NDVI_", RegionFiles[[CombineRun]],
                          ".nc", sep = ""))) {
      print(paste("NDVI raster stack already cropped for: ", RegionFiles[[CombineRun]],
                  sep = ""))
      (next)()
    }
    CombineCDFs(Region = Regions[[CombineRun]], RegionFile = RegionFiles[[CombineRun]],
                Extent = Extents[[CombineRun]])
  } # CombineCDFs
  ##### ERA5 -----
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("HANDLING ERA5 DATA")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  setwd(mainDir)
  source(paste(Dir.Codes, "S2_ERA5.R", sep="/"))
  # kriging ERA5 variable data across study regions for selected time period
  # parallel
  if (Cores > 1) {
    cl <- makeCluster(Cores) # Assuming X node cluster
    registerDoParallel(cl) # Register cores
    for (KrigRegion in 1:length(Regions)) {
      # looping over regions
      print("#################################################")
      print(paste("Kriging ERA5 ", toString(ModVars), " data from across ",
                  RegionFiles[[KrigRegion]], sep = ""))
      foreach(Krigrun = 1:length(ModVars)) %dopar% {
        # looping over variables
        source(paste("Y - Codes", "S0a_Packages.R", sep="/")) # load packages to each core
        source(paste("Y - Codes", "S0b_Directories.R", sep="/")) # set packages for each core
        source(paste(Dir.Codes, "S0c_Functions.R", sep="/")) # Loading miscellaneous functions
        source(paste(Dir.Codes, "S2_ERA5.R", sep="/")) # source function for each core
        ModVars <- c("Tair", "Qsoil1", "Qsoil2", "Qsoil3", "Qsoil4")
        # Checking if this particular data has been kriged already
        if (!file.exists(paste(Dir.ERA.Monthly, "/", ModVars[[Krigrun]],
                               "_mean_", RegionFiles[[KrigRegion]], "_", 1, FromY, "_", 12, To,
                               ".nc", sep = "")))
        {
          sapply(package_vec, install.load.package)
          RasterEra5(Variable = ModVars[[Krigrun]], Region = Regions[[KrigRegion]],
                     RegionFile = RegionFiles[[KrigRegion]], Extent = Extents[[KrigRegion]],
                     FromY = FromY, FromM = 1, ToY = To, ToM = 12, Temporary = "Keep")
        } # check if already kriged
      } # parallel run
    } # Region-loop
    stopCluster(cl) # stop cluster
  } else {
    # non-parallel looping over regions looping over variables Checking if this
    # particular data has been kriged already
    for (KrigRegion in 1:length(Regions)) {
      for (Krigrun in 1:length(ModVars)) {
        if (file.exists(paste(Dir.ERA.Monthly, "/", ModVars[[Krigrun]], "_mean_",
                              RegionFiles[[KrigRegion]], "_", 1, FromY, "_", 12, To, ".nc", sep = ""))) {
          print(paste(ModVars[[Krigrun]], " data already kriged for: ", RegionFiles[[KrigRegion]],
                      sep = ""))
          (next)()
        }
        RasterEra5(Variable = ModVars[[Krigrun]], Region = Regions[[KrigRegion]],
                   RegionFile = RegionFiles[[KrigRegion]], Extent = Extents[[KrigRegion]],
                   FromY = FromY, FromM = 1, ToY = To, ToM = 12, Temporary = "Keep")
      } # Variable-loop
    } # region-loop
  } # RasterEra5 function
  ##### VEGETATION MEMORY -----
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("IDENTIFYING VEGETATION MEMORY")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  setwd(mainDir)
  source(paste(Dir.Codes, "S3_VegetationMemory.R", sep="/"))
  # Calculating vegetation memory
  `%nin%` = Negate(`%in%`) # create a 'not in' statement
  if (Cores > 1) {
    # parallel
    cl <- makeCluster(Cores) # Assuming X node cluster
    registerDoParallel(cl) # registering cores
    for (MemReg in 1:length(Regions)) {
      # looping over regions
      print("#################################################")
      print(paste("Calculating vegetation memory according to ", toString(ModVars),
                  " across ", RegionFiles[[MemReg]], sep = ""))
      foreach(Memrun = 2:length(ModVars)) %dopar% {
        # looping over variables
        source(paste("Y - Codes", "S0a_Packages.R", sep="/")) # load packages to each core
        source(paste("Y - Codes", "S0b_Directories.R", sep="/")) # register directories with each core
        source(paste(Dir.Codes, "S0c_Functions.R", sep="/")) # Loading miscellaneous functions
        source(paste(Dir.Codes, "S3_VegetationMemory.R", sep="/")) # source functions for each core
        ModVars <- c("Tair", "Qsoil1", "Qsoil2", "Qsoil3", "Qsoil4")
        # check if already computed
        if (paste(RegionFiles[[MemReg]], "_Tair_mean-", ModVars[Memrun],
                  "_mean", paste(Lags, collapse = "_"), "_", FromY, "-", To, ".nc",
                  sep = "") %nin% list.files(Dir.Memory)) {
          VegMem(ClimVar = paste(ModVars[Memrun], "_mean", sep = ""), ClimVar2 = "Tair_mean",
                 Region = RegionFiles[[MemReg]], Cumlags = Lags, FromY = FromY,
                 ToY = To)
        }
      } # parallel loop
    } # region-loop
    stopCluster(cl) # stop cluster
  } else {
    # non-parallel looping over regions looping over coefficients of soil layers
    # check if already computed
    for (MemReg in 1:length(Regions)) {
      for (Memrun in 2:length(ModVars)) {
        if (paste(RegionFiles[[MemReg]], "_Tair_mean-", ModVars[Memrun],
                  "_mean", paste(Lags, collapse = "_"), "_", FromY, "-", To, ".nc",
                  sep = "") %nin% list.files(Dir.Memory)) {
          VegMem(ClimVar = paste(ModVars[Memrun], "_mean", sep = ""), ClimVar2 = "Tair_mean",
                 Region = RegionFiles[[MemReg]], Cumlags = Lags, FromY = FromY,
                 ToY = To)
        } else {
          print(paste("Vegetation memory already computed for", ModVars[Memrun],
                      " across:", RegionFiles[[MemReg]], sep = " "))
        }
      } # memory-loop
    } # region-loop
  } # VegMem function
  # scaling coefficients per region to be represented on fixed scales looping over
  # regions
  for (MemReg in 1:length(Regions)) {
    CoeffScaling(ClimVar = ClimVars, ClimVar2 = ClimVars2, Region = list(RegionFiles[[MemReg]],
                                                                         RegionFiles[[MemReg]], 
                                                                         RegionFiles[[MemReg]], 
                                                                         RegionFiles[[MemReg]]),
                 Cumlags = list(Lags, Lags, Lags, Lags), FromY = FromY, ToY = To, UAbs = TRUE)
  } # CoeffScaling function 
} # Fun_Vegetation
####--------------- Fun_PFTs [Regions, RegionFiles, Extents, From, To, Occ]
# (aggregating PFT data, downloading species occurences, building PFT rasters)
# ----
Fun_PFTs <- function(Traits, Regions, RegionFiles, Extents) {
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("OBTAINING DATA FROM BIEN DATA BASE AND CREATING RASTERS")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  source(paste(Dir.Codes, "S4_PFTs.R", sep="/"))
  for(TraitLoop in 1:length(Traits)){
    for(RegionLoop in 1:length(RegionFiles)){
      BIEN(
        Trait = Traits[[TraitLoop]],
        Region = Regions[[RegionLoop]],
        RegionFile = RegionFiles[[RegionLoop]],
        Extent =Extents[[RegionLoop]]
      ) 
    }
  }
} # Fun_PFTs
####--------------- Fun_Compadre [Variables, Regions, RegionFiles, Extents]
# (selecting and preparing data, and making COMPADRE data into rasters) ----
Fun_COMPADRE <- function(Variables, Regions, RegionFiles, Extents) {
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("COMPADRE ANALYSES")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  source(paste(Dir.Codes, "S4_COMPADRE.R", sep="/"))
  # Build rasters of compadre variables across study regions looping over Compadre
  # variables looping over regions
  for (CompVar in 1:length(Variables)) {
    for (CompReg in 1:length(Regions)) {
      if (!file.exists(paste(Dir.Compadre, "/", Variables[[CompVar]], "/",
                             Variables[[CompVar]], "_", RegionFile = RegionFiles[[CompReg]], ".nc",
                             sep = ""))) {
        RasterCOMPADRE(Variable = Variables[[CompVar]], Region = Regions[[CompReg]],
                       RegionFile = RegionFiles[[CompReg]], Extent = Extents[[CompReg]])
      } else {
        print(paste(Variables[[CompVar]], " already rasterised across ",
                    RegionFiles[[CompReg]], sep = ""))
      }
    } # region-loop
  } # CompVar-loop
} # Fun_Compadre
####--------------- FUNCTION CALLS ----
# Fun_Vegetation(Regions = list(c("Portugal", "Spain")),
#                RegionFiles = list("SWEurope"),
#                Extents = list(extent(-10,4.5,35,44)),
#                From = 1982, To = 2015, Lags = 0:12, Cores = 1)

# Fun_PFTs(Traits = list("whole plant height", "leaf nitrogen content per leaf dry mass"),
#          Regions = list("United States of America", "Australia"),
#          RegionFiles = list("Contiguous US", "Australia"),
#          Extents = list(extent(-125, -66, 24, 51), NULL)
# )

Fun_COMPADRE(Variables = list("FastSlow", "Rho", "Pi", 
                              "Reactivity", "FirstStepAtt", 
                              "MaxAmp", "MaxAtt"),
             Regions = list(c("Portugal", "Spain"), "United States of America"),
             RegionFiles = list("SWEurope", "Contiguous US"),
             Extents = list(extent(-10,4.5,35,44), extent(-125, -66, 24, 51)))