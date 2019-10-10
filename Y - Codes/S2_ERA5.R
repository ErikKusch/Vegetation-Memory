###---------------- PREAMBLE ---------------------------------------------------------
Variables_vec <- c("Tair_mean", "Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)",
                        "Soil Moisture (28-100cm)", "Soil Moisture (100-255cm)")
Covariates_vec <- c("Slopes1", "Slopes2", "Slopes3", "Slopes4", "Slopes5", "Slopes6", "Slopes7", "Slopes8",
                    "Slope_aspect_N", "Slope_aspect_E", "Slope_aspect_S", "Slope_aspect_W", "Slope_aspect_U",
                    "Elevation")
Variables_vec <- c("Tair_mean", "Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", 
                        "Soil Moisture (28-100cm)", "Soil Moisture (100-200cm)")
Covariates_vec <- c("Slopes1", "Slopes2", "Slopes3", "Slopes4", "Slopes5", "Slopes6", "Slopes7", "Slopes8", 
                    "Slope_aspect_N", "Slope_aspect_E", "Slope_aspect_S", "Slope_aspect_W", "Slope_aspect_U", 
                    "Elevation")
####--------------- RasterEra5 [Variable, Region, FromY, FromM, ToY, ToM, Temporary]
# (selecting data, downscaling, exporting rasters) ----
RasterEra5 <- function(Variable, Region, RegionFile, Extent, FromY, FromM, ToY, ToM, Temporary){
  VarPos <- which(Variables_vec == paste(Variable, "mean", sep="_")) # position for indexing of variable
  YearVec <- rep(1980:2015, each = 12) # Year vector to indicate months for time frame selection
  # CONSOLE MESSAGE----
  print("#################################################")
  print(paste("Kriging ERA5 ", VariablesNames_vec[VarPos], " data from ", FromM, "/", FromY, " to ", ToM, "/", ToY,
              " across ", RegionFile, sep=""))
  # FORMULAE VECTORS----
  Krig_formula <- "ERA5 ~ Slopes1+Slopes2+Slopes3+Slopes4+Slopes5+Slopes6+Slopes7+Slopes8+Slope_aspect_N+
  Slope_aspect_E+Slope_aspect_S+Slope_aspect_W+Slope_aspect_U+Elevation+
  Slopes1:Slope_aspect_N+Slopes2:Slope_aspect_N+Slopes3:Slope_aspect_N+
  Slopes4:Slope_aspect_N+Slopes5:Slope_aspect_N+Slopes6:Slope_aspect_N+
  Slopes7:Slope_aspect_N+Slopes8:Slope_aspect_N+Slopes1:Slope_aspect_S+
  Slopes2:Slope_aspect_S+Slopes3:Slope_aspect_S+Slopes4:Slope_aspect_S+
  Slopes5:Slope_aspect_S+Slopes6:Slope_aspect_S+Slopes7:Slope_aspect_S+
  Slopes8:Slope_aspect_S"
  # LOAD DATA----
  ## Era5 data
  FirstMonth <- which(YearVec == FromY)[FromM] # first month to consider
  LastMonth <- which(YearVec == ToY)[ToM] # last month to consider
  ras <- list() # create empty list for era5 raster data
  Montquence <- FirstMonth:LastMonth
  ras <- brick(paste(Dir.ERA, "/", Variable ,"_TrainingResolution.nc", sep="")) # loading data
  ras <- ras[[Montquence]] # limitting to sought-after months
  extent(ras) <- c(-180,180,-90,90) # fix extent
  ## Covariates for Kriging
  Cov_coarse <- list() # create empty list
  for(c in 1:length(Covariates_vec)){ # cycle through all covariates and load the data
    Cov_coarse[[c]] <- raster(paste(Dir.KrigCov, "/Co-variates_TrainingResolution.nc", sep=""),
                              varname = Covariates_vec[c])}
  Cov_coarse <- brick(Cov_coarse) # make coarse covariate data into one big brick
  extent(Cov_coarse) <- c(-180,180,-90,90) # fix extent
  Cov_fine <- list() # create empty list
  for(c in 1:length(Covariates_vec)){ # cycle through all covariates and load the data
    Cov_fine[[c]] <- raster(paste(Dir.KrigCov, "/Co-variates_NativeResolution.nc", sep=""),
                            varname = Covariates_vec[c])}
  Cov_fine <- brick(Cov_fine) # make fine covariate data into one big brick
  # REGION SELECTION----
  Shapes <- readOGR(Dir.Mask,'ne_50m_admin_0_countries', verbose = FALSE)
  RegObj <- RegionSelection(Region = Region, RegionFile = RegionFile, Extent = Extent)
  area <- RegObj[[1]]
  location <- RegObj[[2]]
  RegionFile <- RegObj[[3]]
  # CROPPING AND MASKING----
  ## Era5 cropping and masking
  ras <- crop(ras, area) # cropping to extent
  ras <- mask(ras, Shapes[location,]) # masking via Shapefile
  ## Coarse covariate cropping and masking
  Cov_coarse <- crop(Cov_coarse, area) # cropping to extent
  Cov_coarse <- mask(Cov_coarse, Shapes[location,]) # masking via Shapefile
  ## Fine covariate cropping and masking
  Cov_fine <- crop(Cov_fine, area) # cropping to extent
  Cov_fine <- mask(Cov_fine, Shapes[location,]) # masking via Shapefile
  # KRIGING----
  Months1 <- (ToY-FromY-1)*12 # how many months to cover just by years
  Months2 <- abs(ToM-FromM+1) # how many months to cover only taking months of time frame into account
  FullMonths <- Months1+Months2 # total count of months that are covered
  if(FromY == ToY){# if range doesn't exceed a calendar year
    Years <- rep(YearVec[FirstMonth], Months2) }else{
    FromLeft <- 12-FromM+1 # months left in starting year
    ToCovered <- ToM # months to be covered in final year
    Years1 <- rep(FromY, FromLeft)
    Years2 <- rep(ToY, ToCovered)
    if(ToY-FromY > 1){
      Years3 <- rep((FromY+1):(ToY-1), each = 12) # months in full years
      Years <- c(Years1, Years3, Years2) }else{
      Years <- c(Years1, Years2)}}
  Months <- rep(c(1:12), length = length(Years))
  Names <- paste(month.abb, Years, sep="") # combination of month names and years
  ## Preparing Kriging
  Dir.Temp <- paste(Dir.ERA.Monthly,"/Temp_",Variable,"_",RegionFile, sep="")
  dir.create(Dir.Temp)
  TempNames <- paste(rep(YearVec),
                     rep(c("01","02","03","04","05","06","07","08","09","10","11","12")),sep="_")
  # figuring out where to begin with the names
  if(FromM < 10){
    TempStart <- which(TempNames == paste(FromY,"_0",FromM, sep="")) }else{
    TempStart <- which(TempNames == paste(FromY,"_",FromM, sep=""))}
  # figuring out where to stop with the names
  if(ToM < 10){
    TempStop <- which(TempNames == paste(ToY,"_0",ToM, sep="")) }else{
    TempStop <- which(TempNames == paste(ToY,"_",ToM, sep=""))}
  TempNames <- TempNames[TempStart:TempStop]
  Ras_Krig <- list()
  ### Actual Kriging
  counter <- 0
  for(i in 1:length(names(ras))){
    if(paste(TempNames[i], ".nc", sep="") %in% list.files(Dir.Temp)){ # check if this file has already been produced
      print(paste(TempNames[i], "already kriged", sep=" "))
      Ras_Krig[[i]] <- raster(paste(Dir.Temp, "/", TempNames[i], ".nc", sep=""))
      next()}
    counter <- counter + 1
    T_Begin <- Sys.time()
    RasterX <- ras[[i]]# extracting raster from Era5 stack
    # Base and Covariate Coarse Data
    Origin <- as.data.frame(RasterX, xy = TRUE)
    Origin <- na.omit(Origin)
    for(c in 1:length(Covariates_vec)){
      Cov_coarse[[c]][!is.na(RasterX) & is.na(Cov_coarse[[c]])] <- 0 # 0 cells where no info
      Cov_coarse[[c]][is.na(RasterX)] <- NA # ensure same NAs
      Covariate <- as.data.frame(Cov_coarse[[c]], xy = TRUE)
      Covariate <- na.omit(Covariate)
      Origin <- cbind(Origin, Covariate[,3])}
    colnames(Origin) <- c("x","y", "ERA5", Covariates_vec)
    # checking data availability
    for(it_check in 1:length(colnames(Origin))){
      if(length(which(Origin[,it_check] != 0)) < 2){
        stop(paste("The native resolution data does not support kriging using the formula you have specified
                   because ", colnames(Origin)[it_check], " does not contain enough data records for kriging
                   to be performed across the region you have specified (", Region, ").", " You can resolve
                   this issue by either removing the interaction effects containing this variable from the
                   formula or choosing a bigger study region.", sep=""))}}
    OriginK <- Origin
    gridded(OriginK) <- ~x+y
    # CROPPING TARGET
    Cov_fine[[1]][which(is.na(as.vector(Cov_fine[[1]])))] <- 0
    Cov_fine[[1]] <- mask(Cov_fine[[1]], Shapes[location,])
    Target <- as.data.frame(Cov_fine[[1]], xy = TRUE)
    for(c in 2:length(Covariates_vec)){
      Cov_fine[[c]][which(is.na(as.vector(Cov_fine[[c]])))] <- 0
      Cov_fine[[c]] <- mask(Cov_fine[[c]], Shapes[location,])
      Covariate <- as.data.frame(Cov_fine[[c]], xy = TRUE)
      Target <- cbind(Target, Covariate[,3])}
    colnames(Target) <- c("x","y", Covariates_vec)
    TargetK <- Target
    gridded(TargetK) <- ~x+y
    # KRIGING
    invisible(capture.output(
      kriging_result <- autoKrige(
        as.formula(Krig_formula),OriginK, TargetK, verbose = FALSE)))
    Krig_ras <- raster(kriging_result$krige_output)
    Ras_Krig[[i]] <- Krig_ras
    # writing the raster
    writeRaster(Krig_ras, filename = paste(Dir.Temp,"/",TempNames[i],sep=""), overwrite=TRUE, format="CDF")
    if(counter == 1){
      T_End <- Sys.time()
      Duration <- as.numeric(T_End)-as.numeric(T_Begin)
      print(paste("Calculating monthly ERA5 ", VariablesNames_vec[VarPos], " rasters from ", FromM, "/", FromY,
                  " to ", ToM, "/", ToY, " across ", RegionFile, " should finish around: ",
                  as.POSIXlt(T_Begin + Duration*(length(names(ras))-i), tz = Sys.timezone(location=TRUE)), sep=""))
      pb <- txtProgressBar(min = 0, max = length(names(ras)), style = 3)}
    setTxtProgressBar(pb, i)} # kriging loop
  # COMBINING KRIGED ENSEMBLES FROM MEMORY----
  Ras_Krig <- brick(Ras_Krig)
  ras <- Ras_Krig
  # ELIMINATE KRIGING ARTIFACTS OF SOIL MOISTURE BY BOUNDING
  if(Variable == "Qsoil1" | Variable == "Qsoil2" | Variable == "Qsoil3" | Variable == "Qsoil4"){
    values(ras)[which(values(ras) < 0)] <- 0}
  # SAVING DATA----
  setwd(Dir.ERA.Monthly)
  writeRaster(ras, paste(Variable, "_mean_", RegionFile, "_", FromM, FromY, "_", ToM, ToY, sep=""),
              overwrite=TRUE, format="CDF", varname=Variable,
              longname= paste(Variables_vec[VarPos], " mean for years ", FromM, "/",FromY, " to ", ToM,"/",ToY,
                              " across ", Region, sep=""))
  if(Temporary == "Delete"){unlink(Dir.Temp, recursive = TRUE)}
  setwd(mainDir)}# end of RasterEra5 function