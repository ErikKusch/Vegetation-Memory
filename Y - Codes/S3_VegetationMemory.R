Variables_vec <- c("Tair_mean", "Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", 
                        "Soil Moisture (28-100cm)", "Soil Moisture (100-200cm)")
####--------------- VegMem [ClimVar, Region, Cumlags, FromY, ToY] 
# (selecting data, calculating vegetation memory according to specified lags, exporting rasters) ----
VegMem <- function(ClimVar, ClimVar2, Region, Cumlags, FromY, ToY){
  print("#################################################")
  print(paste("Identifying vegetation memory effects of NDVI based on antecedent NDVI and ", 
              VariablesNames_vec[which(Variables_vec == ClimVar2)], ", and ", 
              VariablesNames_vec[which(Variables_vec == ClimVar)], " at lags: ", toString(Cumlags)," across ", 
              Region, sep=""))
  # LOAD DATA----
  ## NDVI/GIMMs
  NDVI_ras <- brick(paste(Dir.Gimms.Monthly,"/NDVI_", Region,".nc", sep=""))
  ## ERA5
  Clim <- list.files(Dir.ERA.Monthly)[grep(pattern = ClimVar, list.files(Dir.ERA.Monthly))] 
  Clim <- Clim[grep(pattern = Region, Clim)] # files in target region
  Clim <- Clim[grep(pattern = FromY, Clim)] # files with correct start date
  Clim <- Clim[grep(pattern = ToY, Clim)] # files with correct end date
  Clim_mean_ras <- brick(paste(Dir.ERA.Monthly, "/", Clim, sep="")) # rasterise
  Clim2 <- list.files(Dir.ERA.Monthly)[grep(pattern = ClimVar2, list.files(Dir.ERA.Monthly))]
  Clim2 <- Clim2[grep(pattern = Region, Clim2)]
  Clim2 <- Clim2[grep(pattern = FromY, Clim2)]
  Clim2 <- Clim2[grep(pattern = ToY, Clim2)]
  Clim2_mean_ras <- brick(paste(Dir.ERA.Monthly, "/", Clim2, sep=""))
  # PREPARE DATA----
  ## Limit NDVI data to ERA5 time frame
  NDVIYears <- rep(1982:2015, each = 12) # Year vector to indicate months for time frame selection
  NDVITo <- max(which(NDVIYears %in% ToY))
  if(min(which(NDVIYears %in% FromY)) == Inf){
    NDVIFrom <- 1 }else{ NDVIFrom <- min(which(NDVIYears %in% FromY))}
  NDVI_ras <- NDVI_ras[[NDVIFrom:NDVITo]] # NDVI data is limitted
  ## Identify data positions
  # establish a mean raster (this sets every cell to NA where any NA is within the time series)
  NATest_ras <- mean(NDVI_ras)
  NATest_vec <- values(NATest_ras) # set values as vector
  Data_Pos <- which(!is.na(NATest_vec)) # select non-NA positions (these are the ones we should build models on)
  # PREPARE RASTERS----
  ModelEval_ras <- NDVI_ras[[1:15]] # select 15 raster layers
  # ModelEval_ras <- NDVI_ras[[1:16]] # select 16 raster layers
  # put names on the layers to tell us what they contain later
  ModelEval_ras <- Fun_NamesRas(raster = ModelEval_ras, ClimVar = ClimVar, ClimVar2 = ClimVar2)
  sink(file = paste(Dir.Memory, "/", Region, ".txt", sep=""))
  print("Working on it")
  sink()
  # MODELS----
  pbi <- 0 # parameter for progress bar
  for(pixel in Data_Pos){ # loop non-NA pixels
    unlink(paste(Dir.Memory, "/", Region, ".txt", sep=""), recursive = TRUE)
    # sink(file = paste(Dir.Memory, "/", Region, ".txt", sep=""))
    # print(paste(pixel, max(Data_Pos), sep="_"))
    sink()
    T_Begin <- Sys.time() # note time when calculation is started (needed for estimation of remaining time)
    ## DATA ----
    ### NDVI stuff -----
    NDVI_vecraw <- as.vector(NDVI_ras[pixel]) # eytract data
    NDVI_vecdet <- detrend(NDVI_vecraw, tt = 'linear') # linear detrending
    # create NDVI data frame
    NDVI_df <- data.frame(Month = rep(1:12, length(NDVI_vecraw)/12), NDVI_raw = NDVI_vecraw, NDVI_de = NDVI_vecdet)
    ## calculate anomalies (Z-scores) and monthly means
    NDVI_df <- transform(NDVI_df, NDVI_Anomalies = ave(NDVI_de, Month, FUN=scale),
                         NDVI_Threshold = ave(NDVI_raw, Month, FUN=function(t) mean(t, na.rm=TRUE)))
    NDVI_df <- NDVI_df[nrow(NDVI_df):1,] # reverse order to read "present to past"
    NDVI_anom <- c(NDVI_df$NDVI_Anomalies, rep(NA, max(Cumlags))) # extract anomalies, adding cumlag NAs
    ThreshPos <- which(NDVI_df$NDVI_Threshold < 0.1) # positions which should be excluded
    if(length(ThreshPos) == length(NDVI_vecraw)){ # if all months should be masked due to NDVImean < 0.1
      # set all in model raster layers to NA for this pixel
      ModelEval_ras[pixel] <- as.numeric(rep(NA, dim(ModelEval_ras)[3])) 
      next()}
    # calculate lag 1
    NDVI_Lag1 <- c(NDVI_anom[-1], NA) # adding one NA for month preceeding data range of NDVI itself
    ### Climate stuff -----
    #### ClimVar ------
    Clim_vec <- as.vector(Clim_mean_ras[pixel]) # extract raw data for pixel (instantenous predictor)
    Clim_vec <- detrend(Clim_vec, tt = 'linear') # linear detrending
    Clim_vec <- Clim_vec[nrow(Clim_vec):1,] # reverse order to read "present to past"
    # calculate cumulative climate indices (antecedent predictor)
    Clim_cum <- rep(NA, length(Cumlags))
    Clim_cum <- as.list(Cumlags)
    position <- 1
    for(lag in Cumlags){
      for(i in 1:(length(Clim_vec)-lag)){
        Clim_cum[[position]] <- c(Clim_cum[[position]], sum(Clim_vec[i:(i+lag)]))}
      Clim_cum[[position]] <- Clim_cum[[position]][-1] # removing initial NA
      # adding enough NAs to bring it up to full length
      Clim_cum[[position]] <- c(Clim_cum[[position]], rep(NA , length(Clim_vec)-length(Clim_cum[[position]])))
      position <- position+1}
    # make data frame of climate stuff
    Clim_df <- as.data.frame(Clim_cum) # make list into data frame
    Clim_df <- cbind(rep(12:1, length(Clim_vec)/12), Clim_vec, Clim_df) # append month index and raw data
    colnames(Clim_df) <- c("Month", "Clim_raw", paste(rep("ClimCum_",length(Cumlags)),Cumlags, sep="")) # column names
    # calculate anomalies
    for(anomaly in 2:length(Clim_df)){# cycle through all the columns of the climate data frame except the month column
      Clim_iter <- with(Clim_df, cbind(Month, Clim_df[,anomaly])) # extract necessary data
      colnames(Clim_iter) <- c("Month", "AnomalyCalc") # set column names
      Clim_iter <- transform(Clim_iter, # calculate anomaly for each month
                             AnomalyCalc = ave(AnomalyCalc, Month, FUN=scale))
      # save to original data frame
      Clim_df[,anomaly] <- Clim_iter$AnomalyCalc}
    #### ClimVar2 ------
    Clim2_vec <- as.vector(Clim2_mean_ras[pixel]) # extract raw data for pixel (instantenous predictor)
    Clim2_vec <- detrend(Clim2_vec, tt = 'linear') # linear detrending
    # calculate cumulative climate indices (antecedent predictor)
    Clim2_cum <- rep(NA, length(Cumlags))
    Clim2_cum <- as.list(Cumlags)
    position <- 1
    for(lag in Cumlags){
      for(i in 1:(length(Clim_vec)-lag)){
        Clim2_cum[[position]] <- c(Clim2_cum[[position]], sum(Clim2_vec[i:(i+lag)]))}
      Clim2_cum[[position]] <- Clim2_cum[[position]][-1] # removing initial NA
      # adding enough NAs to bring it up to full length
      Clim2_cum[[position]] <- c(Clim2_cum[[position]], rep(NA , length(Clim2_vec)-length(Clim2_cum[[position]])))
      position <- position+1}
    # make data frame of climate stuff
    Clim2_df <- as.data.frame(Clim2_cum) # make list into data frame
    Clim2_df <- cbind(rep(12:1, length(Clim2_vec)/12), Clim2_vec, Clim2_df) # append month index and raw data
    colnames(Clim2_df) <- c("Month", "Clim2_raw", paste(rep("Clim2Cum_",length(Cumlags)),Cumlags, sep="")) # column names
    # calculate anomalies
    for(anomaly in 2:length(Clim2_df)){# cycle through all the columns of the climate data frame except the month column
      Clim2_iter <- with(Clim2_df, cbind(Month, Clim2_df[,anomaly])) # extract necessary data
      colnames(Clim2_iter) <- c("Month", "AnomalyCalc") # set column names
      Clim2_iter <- transform(Clim2_iter, # calculate anomaly for each month
                              AnomalyCalc = ave(AnomalyCalc, Month, FUN=scale))
      # save to original data frame
      Clim2_df[,anomaly] <- Clim2_iter$AnomalyCalc}
    ### Combining all the data -----
    ModData_df <- cbind(NDVI_anom, NDVI_Lag1, Clim_df, Clim2_df)
    if(length(ThreshPos) > 0){ # set threshold months to NA if necessary
      ModData_df$NDVI_anom[ThreshPos] <- NA}
    ModData_df <- na.omit(ModData_df) # get rid of NA rows
    if(dim(ModData_df)[1] == 0){ # if there are no anomalies
      ModelEval_ras[pixel] <- c(
        # NA, 
        NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next()
    }
    AICs_df <- data.frame(matrix(rep(NA, length(Cumlags)*length(Cumlags)), nrow = length(Cumlags)))
    counter1 <- 0 # create a counter variable
    ## MODELS ----
    for(ModelIter1 in Cumlags){ # go through all possible cumulative lags
      counter2 <- 0
      for(ModelIter2 in Cumlags){ 
        pca_mat <- as.matrix(cbind(ModData_df$NDVI_Lag1,
                                   ModData_df[, counter1+5], 
                                   ModData_df[, counter2+7+length(Cumlags)]))
        dimnames(pca_mat) <- list(1:length(ModData_df$NDVI_Lag1),                                                                 c("t-1", paste(ClimVar,counter1), paste(ClimVar2, counter2)))
        pca <- rda(pca_mat) # running pca
        ## Extracting PC axes
        pc1 <- summary(pca)[["sites"]][,1]
        pc2 <- summary(pca)[["sites"]][,2]
        pc3 <- summary(pca)[["sites"]][,3]
        Mod <- lm(ModData_df$NDVI_anom ~ pc1 + pc2 + pc3) # full model
        AICs_df[counter1+1, counter2+1] <- AIC(Mod)
        counter2 <- counter2 +1
      }
      counter1 <- counter1 +1
    }
    
    # find best model from data frame  
    for(find in 1:dim(AICs_df)[1]){ # iterating over rows 
      Pos <- which(AICs_df[find,] ==  min(AICs_df))
      if(length(Pos) > 0){
        PCABest <- c(find, Pos) # (clim, clim2)
      }
    }
    
    ### FINAL PCA MODEL ----
    Mod0 <- lm(ModData_df$NDVI_anom ~ 1) # null model
    pca_mat <- as.matrix(cbind(ModData_df$NDVI_Lag1,
                               ModData_df[, PCABest[1]+4], 
                               ModData_df[, PCABest[2]+6+length(Cumlags)]))
    pca <- rda(pca_mat) # running pca
    ## Extracting PC axes
    pc1 <- summary(pca)[["sites"]][,1]
    pc2 <- summary(pca)[["sites"]][,2]
    pc3 <- summary(pca)[["sites"]][,3]
    ## Building models
    Mod0 <- lm(ModData_df$NDVI_anom ~ 1) # null model
    Mod <- lm(ModData_df$NDVI_anom ~ pc1 + pc2 + pc3) # full model
    loadings <- summary(pca)[["species"]] # extract loadings
    coefficients <- Mod$coefficients[2:(dim(summary(pca)[["sites"]])[2]+1)] # extract coefficients
    ## Make coefficients representative by multiplying them with the loadings
    t1newCof <- loadings[1,] * coefficients
    CnewCof <- loadings[2,] * coefficients
    C2newCof <- loadings[3,] * coefficients
    ## Saving information to vectors
    coeffst1 <- sum(t1newCof) # NDVI-1
    coeffsC <- sum(CnewCof) # ClimVar
    coeffsC2 <- sum(C2newCof) # ClimVar2
    c_NDVI <- coeffst1 # ndvi coefficient
    c_Clim <- coeffsC # climate coefficient
    c_Clim2 <- coeffsC2 # climate 2 coefficient
    if(anova(Mod0, Mod)$RSS[1] > anova(Mod0, Mod)$RSS[2]){ # only save p value if model is an improvement
      ps <- anova(Mod0, Mod)$'Pr(>F)'[2] 
    }else{ # if model is not an improvement over null, set p to 1
      ps <- 1}
    ## EXPLAINED VARIANCE----
    VarPart <- with(ModData_df, modEvA::varPart(A = summary(lm(NDVI_anom ~ NDVI_Lag1))[["r.squared"]], 
                                                B = summary(lm(NDVI_anom ~ ModData_df[, PCABest[1]+4]))[["r.squared"]], 
                                                C = summary(lm(NDVI_anom ~ ModData_df[, PCABest[2]+6+length(Cumlags)]))[["r.squared"]], 
                                                AB = summary(lm(NDVI_anom ~ NDVI_Lag1+ModData_df[, PCABest[1]+4]))[["r.squared"]],
                                                BC = summary(lm(NDVI_anom ~ ModData_df[, PCABest[1]+4]))[["r.squared"]], 
                                                AC = summary(lm(NDVI_anom ~ NDVI_Lag1+ModData_df[, PCABest[2]+6+length(Cumlags)]))[["r.squared"]], 
                                                ABC = summary(lm(NDVI_anom ~ NDVI_Lag1+ModData_df[, PCABest[1]+4]+ModData_df[, PCABest[2]+6+length(Cumlags)]))[["r.squared"]], 
                                                A.name = "t-1", B.name = "Qsoil1", C.name = "Tair", main = "Memory Components", plot = FALSE)
    )
    Vars <- c(1-VarPart[8,1], # total variance
              VarPart[1,1], # t-1
              VarPart[2,1], # Qsoil1 
              VarPart[3,1], # Tair
              VarPart[4,1], # t-1 + qsoil1
              VarPart[6,1], # t-1 + Tair
              VarPart[5,1], # Qsoil1 + Tair
              VarPart[7,1]) # shared by all
    Vars[which(Vars < 0)] <- 0
    
    # with(ModData_df, varpart(Y = ModData_df, X = NDVI_anom ~ NDVI_Lag1, ~ModData_df[, PCABest[1]+4], ~ModData_df[, PCABest[2]+6+length(Cumlags)], data = ModData_df))
    
    ## WRITING INFORMATION TO RASTERS----
    ModelEval_ras[pixel] <- as.numeric(c(
      # summary(Mod)[["adj.r.squared"]], 
      AIC(Mod), c_NDVI, c_Clim, PCABest[1]-1, c_Clim2, PCABest[2]-1, 
                                         Vars, ps)) # saving model information to raster
    ## Updating progress bar----
    if(exists("pb") == FALSE){ # if we are currently on the first pixel
      T_End <- Sys.time() # note end time
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # calculate the time it took to establish and select models
      ## Put an estimator up on the console that tells the user when to expect the program to finish its current run
      print(paste("Calculating Vegetation Memory effects across ", Region, " should take ", Duration*length(Data_Pos)/3600, " and finish around: ", 
                  as.POSIXlt(T_Begin + Duration*length(Data_Pos), tz = Sys.timezone(location=FALSE)), sep=""))
      ## Update progress bar
      pb <- txtProgressBar(min = 0, max = length(Data_Pos), style = 3)
    }
    pbi <- pbi + 1 ## Update progress bar
    setTxtProgressBar(pb, pbi)} # end of pixel loop
  unlink(paste(Dir.Memory, "/", Region, ".txt", sep=""), recursive = TRUE)
  ### Save data ----- 
  writeRaster(ModelEval_ras, filename = paste(Dir.Memory,"/", Region, "_", ClimVar2, "-", ClimVar, 
                                              paste(Cumlags, collapse="_"),"_",FromY,"-", ToY, ".nc",sep=""),
              overwrite=TRUE, format="CDF")
  setwd(mainDir)
}# end of VegMem function 