###---------------- PREAMBLE ---------------------------------------------------------
Variables_vec <- c("Tair_mean", "Qsoil1_mean")
Download_vec <- c("2m_temperature", "volumetric_soil_water_layer_1")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)")
# Covariates_vec <- c("Slopes1", "Slopes2", "Slopes3", "Slopes4", "Slopes5", "Slopes6", "Slopes7", "Slopes8",
#                     "Slope_aspect_N", "Slope_aspect_E", "Slope_aspect_S", "Slope_aspect_W", "Slope_aspect_U",
#                     "Elevation")
####--------------- RasterEra5 [Variable, Region, FromY, FromM, ToY, ToM, Temporary]
# (selecting data, downscaling, exporting rasters) ----
RasterEra5 <- function(Variable, Region, RegionFile, FromY, FromM, ToY, To, ToM){
  VarPos <- which(Variables_vec == paste(Variable, "mean", sep="_")) # position for indexing of variable
  YearVec <- rep(FromY:ToY, each = 12) # Year vector to indicate months for time frame selection
  Dir.eraiter <- file.path(Dir.ERA, RegionFile)
  dir.create(Dir.eraiter)
  # CONSOLE MESSAGE----
  print("#################################################")
  print(paste("Obtaining ERA5-Land ", VariablesNames_vec[VarPos], " data from ", FromM, "/", FromY, " to ", ToM, "/", ToY,
              " across ", RegionFile, sep=""))
  
  output <- KrigR::download_ERA(Variable = Download_vec[VarPos],
                      DateStart = paste(FromY, str_pad(FromM, 2, side = "left", 0), "01", sep = "-"),
                      DateStop = paste(ToY, str_pad(ToM, 2, side = "left", 0), "31", sep = "-"),
                      Extent = Region,
                      Dir = Dir.eraiter,
                      FileName = paste(Variable, RegionFile, sep = "_"),
                      API_User = API_User,
                      API_Key = API_Key
                      )
  return(output)
  }# end of RasterEra5 function