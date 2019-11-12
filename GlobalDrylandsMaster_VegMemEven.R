rm(list = ls()) # clearing environment
####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions
####--------------- VARIABLE VECTORS ----
ModVars <- c("Tair", "Qsoil1")
             # , "Qsoil2", "Qsoil3", "Qsoil4")
ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
###---------------- FUNCTIONS ---------------------------------------------------------
####--------------- Fun_Vegetation [Regions, RegionFiles, Extents, From, To, Lags, Cores]
# (selecting and preparing data, and calculating vegetation memory) ----
Fun_Vegetation <- function(Regions, RegionFiles, Extents, From, To, Lags, Cores) {
  FromY <- (From - ceiling(1/12 * max(Lags))) # Figuring out real start year after factoring in lags
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
} # Fun_Vegetation

####--------------- FUNCTION CALLS ----
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
Extents <- Extents[seq(2,Tiles,2)]

Fun_Vegetation(Regions = as.list(rep("Drylands", length(Extents))),
               RegionFiles = as.list(paste("Drylands_", seq(2,Tiles,2), sep="")),
               Extents = Extents,
               From = 1982, To = 2015, Lags = 0:12, Cores = 2)