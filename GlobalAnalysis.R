rm(list = ls()) # clearing environment
####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions
####--------------- APISettings ----
if(file.exists("Y - Codes/PersonalSettings.R")){
  source("Y - Codes/PersonalSettings.R") # Loading miscellaneous functions 
}else{
  API_User <- readline("ECMWF API User:")
  API_Key <- readline("ECMWF API Key:")
}
####--------------- VARIABLE VECTORS ----
ModVars <- c("Tair", "Qsoil1")
# ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
# ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
###---------------- FUNCTIONS ---------------------------------------------------------
####--------------- GlobalDrylands [numberofCores]
# (Breaking up global dryland shapefile into bands and enact vegetation memory function) ----
GlobalDrylands <- function(numberOfCores = 13){
  cl <- makeCluster(numberOfCores) # Assuming X node cluster
  registerDoParallel(cl) # registering cores
  setwd(Dir.Mask)
  Drylands <- shapefile("dryland_2")
  Extents <- list()
  Tiles <- floor((abs(extent(Drylands)[4])+abs(extent(Drylands)[3]))/10)
  for(i in 0:Tiles){
    Extents[[i+1]] <- extent(extent(Drylands)[1], extent(Drylands)[2],
                             extent(Drylands)[4]-10*(i+1),
                             extent(Drylands)[4]-10*i)
  }
  setwd(mainDir)
  Regions <- as.list(rep("Drylands", length(Extents)))
  RegionFiles <- as.list(paste("Drylands_", seq(1,Tiles+1,1), sep=""))
  
  foreach(Para = 1:length(Regions)) %dopar% {
    source(paste("Y - Codes", "S0a_Packages.R", sep="/")) # load packages to each core
    source(paste("Y - Codes", "S0b_Directories.R", sep="/")) # set packages for each core
    source(paste(Dir.Codes, "S0c_Functions.R", sep="/")) # Loading misc functions
    if(file.exists("Y - Codes/PersonalSettings.R")){
      source("Y - Codes/PersonalSettings.R") # Loading miscellaneous functions 
    }else{
      API_User <- readline("ECMWF API User:")
      API_Key <- readline("ECMWF API Key:")
    }
    
    ####--------------- Fun_Vegetation [Regions, RegionFiles, Extents, From, To, Lags, Cores]
    # (selecting and preparing data, and calculating vegetation memory) ----
    Fun_Vegetation <- function(Regions, RegionFiles, Extents, From, To, Lags, Cores, Para) {
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
      Region_ras <- raster(paste(Dir.Gimms.Monthly, "/NDVI_", RegionFiles, ".nc", sep = ""))
      
      ##### ERA5 -----
      print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      print("HANDLING ERA5 DATA")
      print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
      setwd(mainDir)
      source(paste(Dir.Codes, "S2_ERA5.R", sep="/"))
      
      for(Krigrun in 1:length(ModVars)){
        if(!file.exists(file.path(Dir.ERA, paste0(ModVars[Krigrun], "_", RegionFiles, ".nc")))){
          ResampERA <- RasterEra5(Variable = ModVars[[Krigrun]], 
                                  Region = extent(Region_ras),
                                  RegionFile = RegionFiles,
                                  FromY = FromY, 
                                  FromM = 1, 
                                  ToY = To, 
                                  ToM = 12) 
          ResampERA2 <- resample(ResampERA, Region_ras) 
          ResampERA3 <- mask(ResampERA2, crop(Drylands, extent(Region_ras)))
          writeRaster(x = ResampERA3, filename = file.path(Dir.ERA, paste0(ModVars[Krigrun], "_", RegionFiles, ".nc")),
                      format = "CDF", overwrite = TRUE)
        }else{
          print(paste("ERA5-Land data already downloaded for", ModVars[Krigrun]))
        }
      }
      try(unlink(file.path(Dir.ERA, RegionFiles), recursive = TRUE))
    
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
    } # Fun_Vegetation

    ModVars <- c("Tair", "Qsoil1")
    # ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
    # ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
    Fun_Vegetation(Regions = Regions[[Para]],
                   RegionFiles = RegionFiles[[Para]],
                   Extents = list(Extents[[Para]]),
                   From = 1982, To = 2015, Lags = 0:12, Cores = 1, Para = Para)
  }
  stopCluster(cl)
}

####--------------- Fun_Plots [Variables, Regions, RegionFiles, Legends, SoilLayers]
# (selecting and preparing data, and making COMPADRE data into rasters) ----
Fun_Plots <- function(Region, Scaled){
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("PRODUCING PLOTS OF VEGETATION MEMORY")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  source(paste(Dir.Codes, "SX_Plots.R", sep="/"))
  Fun_Plot(Region = Region, Scaled = Scaled)
} # Fun_Plots

####--------------- FUNCTION CALLS ----
GlobalDrylands(numberOfCores = 4
  # numberOfCores = detectCores()
  )

## making a global raster of memory effects
setwd(Dir.Memory)
files <- list.files()
files <- files[grep(files, pattern = "2015.nc")][-1] # not using Drylands_1 because of data errors
ls <- list()# loading data
for(i in 1:length(files)){
  ls[[i]] <- brick(files[i])
}
Big_ras <- ls[[1]] # merging data
for(i in 2:length(ls)){
  Big_ras <- merge(Big_ras, ls[[i]])
}
writeRaster(Big_ras, "GlobalDrylands.nc", format="CDF") # saving data
rm(ls)
rm(Big_ras)

Fun_Plots(Region = "GlobalDrylands", Scaled = TRUE)