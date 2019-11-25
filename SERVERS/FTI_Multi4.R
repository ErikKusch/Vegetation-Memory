rm(list = ls()) # clearing environment
####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions
####--------------- VARIABLE VECTORS ----
ModVars <- c("Tair", "Qsoil1")
ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
###---------------- FUNCTIONS ---------------------------------------------------------
numberOfCores = 4
cl <- makeCluster(numberOfCores) # Assuming X node cluster
registerDoParallel(cl) # registering cores

foreach(Para = 13:16) %dopar% {
  outputFile <-file(paste("output", Para,".txt", sep=""))
  
  library(raster)
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
  
  Regions = as.list(rep("Drylands", length(Extents)))
  RegionFiles = as.list(paste("Drylands_", seq(1,Tiles+1,1), sep=""))
  Extents = Extents
  
  Indicators <- rep(c(2:5), each = 4)
  froms <- rep(c(1982, 1993, 2001, 2009), 4)
  tos <- rep(c(1991, 1999, 2007, 2015), 4)
  
  ModVars <- c("Tair", "Qsoil1")
  ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
  ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")
  Lags <- 0:12
  FromY <- (froms[Para] - ceiling(1/12 * max(Lags))) # Figuring out real start year after factoring in lags
  ##### ERA5 -----
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("HANDLING ERA5 DATA")
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  setwd(mainDir)
  source(paste("Y - Codes", "S0a_Packages.R", sep="/")) # load packages to each core
  source(paste("Y - Codes", "S0b_Directories.R", sep="/")) # set packages for each core
  source(paste(Dir.Codes, "S0c_Functions.R", sep="/")) # Loading
  source(paste(Dir.Codes, "S2_ERA5.R", sep="/"))
  # QSOIL
  Krigrun <- 2
  if (file.exists(paste(Dir.ERA.Monthly, "/", ModVars[[Krigrun]], "_mean_",
                        RegionFiles[[Indicators[Para]]], "_", 1, FromY, "_", 12, tos[Para], ".nc", sep = ""))) {
    print(paste(ModVars[[Krigrun]], " data already kriged for: ", RegionFiles[[Indicators[Para]]],
                sep = ""))
  }else{
    tryCatch({
      RasterEra5(Variable = ModVars[[Krigrun]], 
                 Region = Regions[[Indicators[Para]]],
                 RegionFile = RegionFiles[[Indicators[Para]]], 
                 Extent = Extents[[Indicators[Para]]],
                 FromY = FromY, FromM = 1, 
                 ToY = tos[Para], ToM = 12, 
                 Temporary = "Keep")
    }, error = function(e) {
      writeLines(Para, outputFile)
    })
    -----------------------------
      close(outputFile)
  }
}
stopCluster(cl)