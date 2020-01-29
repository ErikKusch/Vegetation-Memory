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
Variables_vec <- c("Tair_mean", "Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", 
                        "Soil Moisture (28-100cm)", "Soil Moisture (100-200cm)")

##### ----------------- DEFINING PARAMETERS --------------
ClimVar <- "Tair_mean"
ClimVar2 <- "Qsoil1_mean"
Region <- "GlobalDrylands"
Cumlags <- 0:12
FromY <- 1981
ToY <- 2015

##### ----------------- CREATING GLOBAL RASTERS OF DATA --------------
`%nin%` = Negate(`%in%`) # create a 'not in' statement


###### ---- NDVI ------
setwd(Dir.Gimms.Monthly)
NCs <- list.files(pattern = "Drylands_")

if("GlobalDrylandsNDVI.nc" %nin% list.files()){
  NDVIRasters <- list()
  for(Band in 1:length(NCs)){
    NDVIRasters[[Band]] <- brick(NCs[Band])
  }
  Merged <- merge(NDVIRasters[[1]], NDVIRasters[[2]], NDVIRasters[[3]], NDVIRasters[[4]], 
                  NDVIRasters[[5]], NDVIRasters[[6]], NDVIRasters[[7]], NDVIRasters[[8]], 
                  NDVIRasters[[9]], NDVIRasters[[10]], NDVIRasters[[11]], NDVIRasters[[12]], 
                  NDVIRasters[[13]])
  writeRaster(Merged, "GlobalDrylandsNDVI",
              overwrite=TRUE, format="CDF", varname="NDVI",
              longname= paste("Monthly NDVI mean for years ", "01", "/", FromY+1, 
                              " to ", "12", "/", ToY, " across ", "Global Drylands", sep="")) 
}else{
  print("Loading NDVI file")
  Merged <- brick("GlobalDrylandsNDVI.nc")
}

NDVI_ras <- Merged

###### ---- ABIOTIC ------
setwd(Dir.ERA.Monthly)
NCs <- list.files(pattern = ".nc")

GlobalRasters <- list()
for(Var in 1:length(ModVars)){
  VarRasters <- list()
  Bands <- NCs[grep(NCs, pattern = ModVars[Var])]
  if(paste(ModVars[[Var]], "_mean_", "GlobalDrylands", "_", 
           "01", FromY, "_", "12", ToY, ".nc", sep="") %nin% list.files()){ # if global drylands variable raster isn't there yet
    for(Band in 1:length(Bands)){
      VarRasters[[Band]] <- brick(Bands[Band])
    }
    Merged <- merge(VarRasters[[1]], VarRasters[[2]], VarRasters[[3]], VarRasters[[4]], 
                    VarRasters[[5]], VarRasters[[6]], VarRasters[[7]], VarRasters[[8]], 
                    VarRasters[[9]], VarRasters[[10]], VarRasters[[11]], VarRasters[[12]], 
                    VarRasters[[13]])
    writeRaster(Merged, paste(ModVars[[Var]], "_mean_", "GlobalDrylands", "_", 
                              "01", FromY, "_", "12", ToY, sep=""),
                overwrite=TRUE, format="CDF", varname=ModVars[[Var]],
                longname= paste(ModVars[[Var]], " mean for years ", "01", "/", FromY, 
                                " to ", "12", "/", ToY, " across ", "Global Drylands", sep="")) 
  }else{
    print(paste("Loading", ModVars[Var] ,"file"))
    Merged <- brick(paste(ModVars[[Var]], "_mean_", "GlobalDrylands", "_", 
                          "01", FromY, "_", "12", "2015.nc", sep=""))
  }
  GlobalRasters[[Var]] <- Merged
}

Clim_mean_ras <- GlobalRasters[[1]]
Clim2_mean_ras <- GlobalRasters[[2]]

##### ----------------- PREPARING MODELS --------------
## identifying a random set of global cells which to analyse ----
print("Calculating possible analysis cells. This might take a while.")
MeanAir <- mean(Clim_mean_ras)
ModelCells <- which(!is.na(values(MeanAir)))
set.seed(42)
ModCells <- sample(x = ModelCells, size = 10000)

## Data Preparation ----
Compare_df <- data.frame( # data frame for final comparison
  AIC = NA, # AIC of final model, comes up 3 times
  Lag = NA, # comes up 2 times + eone time as 1 for AR1
  Effect = NA,  # 3 values
  Variable = NA, # one value for each variable
  Method = NA, # how lags were selected
  Cell = NA # which cell we analysed
)

##### !!!!!!! MODELS !!!!!! ---------------------------------------------------------
for(pixel in ModCells){
  if(pixel == ModCells[1]){
    T_Begin <- Sys.time() # note time when calculation is started (needed for estimation of remaining time) 
  }else{
    pbi <- pbi + 1 ## Update progress bar
    setTxtProgressBar(pb, pbi)
  }
  ##### ----------------- LAG DATA PREPARATION --------------
  ### NDVI stuff -----
  NDVI_vecraw <- as.vector(NDVI_ras[pixel]) # eytract data
  if(is.na(sum(NDVI_vecraw))){
    next()
  }
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
    # print("Threshold Criteria triggered")
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
  
  ##### ----------------- MODELS LINEAR-SINGLE SELECTION --------------
  Mods_ls <- as.list(rep(NA, length(Cumlags))) # List of models
  for(ModelIter in Cumlags){ # go through all possible cumulative lags
    Mod <- lm(NDVI_anom ~ Clim_df[,ModelIter+3]) # full model
    Mods_ls[[ModelIter+1]] <- Mod # save model to list of models
  }
  ### Selecting best model -----
  AICs <- sapply(X = Mods_ls, FUN = AIC) # calculate AICs for each model
  BestC1 <- which(abs(AICs) == min(abs(AICs), na.rm = TRUE))[1] # best model, if same value present use first
  Mods_ls2 <- as.list(rep(NA, length(Cumlags))) # List of models
  for(ModelIter in Cumlags){ # go through all possible cumulative lags
    Mod <- lm(NDVI_anom ~ Clim2_df[,ModelIter+3]) # full model
    Mods_ls2[[ModelIter+1]] <- Mod # save model to list of models
  }
  ### Selecting best model -----
  AICs2 <- sapply(X = Mods_ls2, FUN = AIC) # calculate AICs for each model
  BestC2 <- which(abs(AICs2) == min(abs(AICs2), na.rm = TRUE))[1] # best model, if same value present use first
  ### Combining all the data -----
  ModData_df <- data.frame(NDVI_anom = NDVI_anom, 
                           NDVI_Lag1 = NDVI_Lag1, 
                           C1 = Clim_df[,BestC1+2], 
                           C2 = Clim2_df[,BestC2+2])
  if(length(ThreshPos) > 0){ # set threshold months to NA if necessary
    ModData_df$NDVI_anom[ThreshPos] <- NA}
  ModData_df <- na.omit(ModData_df) # get rid of NA rows
  if(dim(ModData_df)[1] == 0){ # if there are no anomalies
    # print("Threshold Criteria triggered")
    next()
  }
  ### Establishing models----
  ## PCA of our variables
  pca_mat <- matrix(cbind(ModData_df$NDVI_Lag1,ModData_df$C1, ModData_df$C2), 
                    ncol = 3, byrow = FALSE, dimnames = list(1:length(ModData_df$NDVI_Lag1), 
                                                             c("t-1", ClimVar, ClimVar2))) # pca matrix
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
  if(anova(Mod0, Mod)$RSS[1] > anova(Mod0, Mod)$RSS[2]){ # only save p value if model is an imporvement
    ps <- anova(Mod0, Mod)$'Pr(>F)'[2] 
  }else{ # if model is not an improvement over null, set p to 1
    ps <- 1}
  ### Selecting best model -----
  # Best <- which(abs(coeffsC) == max(abs(coeffsC)))[1] # select fro strongest soil memory effect
  c_NDVI <- coeffst1 # ndvi coefficient
  c_Clim <- coeffsC # climate coefficient
  c_Clim2 <- coeffsC2 # climate 2 coefficient
  p_Mod <- ps # p-value is set to p-value of best model
  AICMod <- AIC(Mod) # AIC is set to best AIC
  BestlagC1 <- BestC1-1  # this is the lag at which best model was observed
  BestlagC2 <- BestC2-1  # this is the lag at which best model was observed
  
  ### Saving info to plotting frame -----
  Compare_dfBIND <- data.frame( # data frame for final comparison
    AIC = c(AIC(Mod), AIC(Mod), AIC(Mod)), # AIC of final model, comes up 3 times
    Lag = c(1, BestlagC1, BestlagC2), # comes up 2 times + eone time as 1 for AR1
    Effect = c(c_NDVI, c_Clim, c_Clim2),  # 3 values
    Variable = c("t-1", "Qsoil1", "Tair"), # one value for each variable
    Method = rep("Linear-Single", 3), # how lags were selected
    Cell = rep(pixel, 3) # which cell we analysed
  )
  Compare_df <- rbind(Compare_df, Compare_dfBIND)
  
  ##### ----------------- MODELS LINEAR-COMBINED SELECTION --------------
  Mods_ls <- as.list(rep(NA, length(Cumlags))) # List of models
  for(ModelIter in Cumlags){ # go through all possible cumulative lags
    Mod <- lm(NDVI_anom ~ NDVI_Lag1 + Clim_df[,ModelIter+3]) # full model
    Mods_ls[[ModelIter+1]] <- Mod # save model to list of models
  }
  ### Selecting best model -----
  AICs <- sapply(X = Mods_ls, FUN = AIC) # calculate AICs for each model
  BestC1 <- which(abs(AICs) == min(abs(AICs), na.rm = TRUE))[1] # best model, if same value present use first
  Mods_ls2 <- as.list(rep(NA, length(Cumlags))) # List of models
  for(ModelIter in Cumlags){ # go through all possible cumulative lags
    Mod <- lm(NDVI_anom ~ NDVI_Lag1 + Clim2_df[,ModelIter+3]) # full model
    Mods_ls2[[ModelIter+1]] <- Mod # save model to list of models
  }
  ### Selecting best model -----
  AICs2 <- sapply(X = Mods_ls2, FUN = AIC) # calculate AICs for each model
  BestC2 <- which(abs(AICs2) == min(abs(AICs2), na.rm = TRUE))[1] # best model, if same value present use first
  ### Combining all the data -----
  ModData_df <- data.frame(NDVI_anom = NDVI_anom, 
                           NDVI_Lag1 = NDVI_Lag1, 
                           C1 = Clim_df[,BestC1+2], 
                           C2 = Clim2_df[,BestC2+2])
  if(length(ThreshPos) > 0){ # set threshold months to NA if necessary
    ModData_df$NDVI_anom[ThreshPos] <- NA}
  ModData_df <- na.omit(ModData_df) # get rid of NA rows
  if(dim(ModData_df)[1] == 0){ # if there are no anomalies
    # print("Threshold Criteria triggered")
    next()
  }
  ### Establishing models----
  ## PCA of our variables
  pca_mat <- matrix(cbind(ModData_df$NDVI_Lag1,ModData_df$C1, ModData_df$C2), 
                    ncol = 3, byrow = FALSE, dimnames = list(1:length(ModData_df$NDVI_Lag1), 
                                                             c("t-1", ClimVar, ClimVar2))) # pca matrix
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
  if(anova(Mod0, Mod)$RSS[1] > anova(Mod0, Mod)$RSS[2]){ # only save p value if model is an imporvement
    ps <- anova(Mod0, Mod)$'Pr(>F)'[2] 
  }else{ # if model is not an improvement over null, set p to 1
    ps <- 1}
  ### Selecting best model -----
  # Best <- which(abs(coeffsC) == max(abs(coeffsC)))[1] # select fro strongest soil memory effect
  c_NDVI <- coeffst1 # ndvi coefficient
  c_Clim <- coeffsC # climate coefficient
  c_Clim2 <- coeffsC2 # climate 2 coefficient
  p_Mod <- ps # p-value is set to p-value of best model
  AICMod <- AIC(Mod) # AIC is set to best AIC
  BestlagC1 <- BestC1-1  # this is the lag at which best model was observed
  BestlagC2 <- BestC2-1  # this is the lag at which best model was observed
  
  ### Saving info to plotting frame -----
  Compare_dfBIND <- data.frame( # data frame for final comparison
    AIC = c(AIC(Mod), AIC(Mod), AIC(Mod)), # AIC of final model, comes up 3 times
    Lag = c(1, BestlagC1, BestlagC2), # comes up 2 times + eone time as 1 for AR1
    Effect = c(c_NDVI, c_Clim, c_Clim2),  # 3 values
    Variable = c("t-1", "Qsoil1", "Tair"), # one value for each variable
    Method = rep("Linear-Combined", 3), # how lags were selected
    Cell = rep(pixel, 3) # which cell we analysed
  )
  Compare_df <- rbind(Compare_df, Compare_dfBIND)
  
  ##### ----------------- MODELS FULLY COMBINED PCA SELECTION --------------
  ModData_df <- cbind(NDVI_anom, NDVI_Lag1, Clim_df, Clim2_df)
  if(length(ThreshPos) > 0){ # set threshold months to NA if necessary
    ModData_df$NDVI_anom[ThreshPos] <- NA}
  ModData_df <- na.omit(ModData_df) # get rid of NA rows
  if(dim(ModData_df)[1] == 0){ # if there are no anomalies
    # print("Threshold Criteria triggered")
    next()
  }
  AICs_df <- data.frame(matrix(rep(NA, length(Cumlags)*length(Cumlags)), nrow = length(Cumlags)))
  counter1 <- 0 # create a counter variable
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
  
  ### FINAL PCA MODEL
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
  
  ### Saving info to plotting frame -----
  Compare_dfBIND <- data.frame( # data frame for final comparison
    AIC = c(AIC(Mod), AIC(Mod), AIC(Mod)), # AIC of final model, comes up 3 times
    Lag = c(1, PCABest[1], PCABest[2]), # comes up 2 times + eone time as 1 for AR1
    Effect = c(c_NDVI, c_Clim, c_Clim2),  # 3 values
    Variable = c("t-1", "Qsoil1", "Tair"), # one value for each variable
    Method = rep("PCA-Combined", 3), # how lags were selected
    Cell = rep(pixel, 3) # which cell we analysed
  )
  Compare_df <- rbind(Compare_df, Compare_dfBIND)
  
  
  
  ## Updating progress bar----
  if(pixel == ModCells[1]){ # if we are currently on the first pixel
    T_End <- Sys.time() # note end time
    Duration <- as.numeric(T_End)-as.numeric(T_Begin) # calculate the time it took to establish and select models
    ## Put an estimator up on the console that tells the user when to expect the program to finish its current run
    print(paste("Calculating Vegetation Memory effects across ", Region, " should finish around: ", 
                as.POSIXlt(T_Begin + Duration*length(ModCells), tz = Sys.timezone(location=FALSE)), sep=""))
  ## Update progress bar
  pb <- txtProgressBar(min = 0, max = length(ModCells), style = 3)
  pbi <- 0
  pbi <- pbi + 1 ## Update progress bar
  setTxtProgressBar(pb, pbi)
  }
} # !!!!!! MODELS !!!!!!


##### PLOTTING COMPARISONS ---------------------------------------------------------
Compare_df <- na.omit(Compare_df)
ggplot(aes(x = Variable, y = AIC, col = Method), data = Compare_df) + 
  geom_boxplot() + theme_bw() + ggtitle("Feature Selection")

# Plot_df <- Compare_df
# Plot_df <- Plot_df[which(Compare_df$Variable == "Tair"), ]
# Plot_df <- Plot_df[with(Plot_df, order(Cell)), ]
# Plot_df$Cells <- sort(Plot_df$Cell)
# 
# ggplot(aes(x = Cells, y = AIC, col = Method), data = Plot_df) + 
#   geom_line() + theme_bw() + ggtitle("Tair Feature Selection")
setwd(mainDir)
