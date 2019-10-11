Variables_vec <- c("Tair_mean", "Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
VariablesNames_vec <- c("Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", 
                        "Soil Moisture (28-100cm)", "Soil Moisture (100-200cm)")
####--------------- VegMem [ClimVar, Region, Cumlags, FromY, ToY] 
# (selecting data, calculating vegetation memory according to specified lags, exporting rasters) ----
VegMem <- function(ClimVar, ClimVar2, Region, Cumlags, FromY, ToY){
  print("#################################################")
  print(paste("Identifying vegetation memory effects of NDVI based on antecedent NDVI and ", 
              VariablesNames_vec[which(Variables_vec == ClimVar2)], " (immediate effects) and ", 
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
  ModelEval_ras <- NDVI_ras[[1:10]] # select six raster layers
  # put names on the layers to tell us what they contain later
  ModelEval_ras <- Fun_NamesRas(raster = ModelEval_ras, ClimVar = ClimVar, ClimVar2 = ClimVar2)
  # MODELS----
  for(pixel in Data_Pos){ # loop non-NA pixels
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
    Clim2_df <- data.frame(Month = rep(1:12, length(Clim2_vec)/12), Clim2_raw = Clim2_vec)
    # calculate anomalies
    Clim2_df <- transform(Clim2_df,
                          Clim2_Anomalies = ave(Clim2_raw, Month, FUN=scale))
    Clim2_df <- Clim2_df[nrow(Clim2_df):1,] # reverse order to read "present to past"
    Clim2_anom <- Clim2_df$Clim2_Anomalies
    ### Combining all the data -----
    ModData_df <- cbind(NDVI_anom[1:length(Clim2_anom)], NDVI_Lag1[1:length(Clim2_anom)], Clim_df, Clim2_anom)
    if(length(ThreshPos) > 0){ # set threshold months to NA if necessary
      ModData_df$NDVI_anom[ThreshPos] <- NA}
    ModData_df <- na.omit(ModData_df) # get rid of NA rows
    if(dim(ModData_df)[1] == 0){ # if there are no anomalies
      ModelEval_ras[pixel] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      next()
    }
    ## MODELS ----
    ### Establishing models----
    # list to save Model objects
    Mods_ls <- as.list(rep(NA, length(Cumlags))) # List of models
    ps <- rep(NA, length(Cumlags)) # p-values
    coeffst1 <- rep(NA, length(Cumlags)) # coefficients of NDVI-1
    coeffsC <- rep(NA, length(Cumlags)) # coefficients of ClimVar
    coeffsC2 <- rep(NA, length(Cumlags)) # coefficients of ClimVar2
    # itterate over all climate lags
    counter <- 0 # create a counter variable
    for(ModelIter in Cumlags){ # go through all possible cumulative lags
      ## PCA of our variables
      pca_mat <- matrix(cbind(ModData_df$NDVI_Lag1,ModData_df[, counter+4], ModData_df$Clim2_anom), 
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
      coeffst1[counter+1] <- sum(t1newCof) # NDVI-1
      coeffsC[counter + 1] <- sum(CnewCof) # ClimVar
      coeffsC2[counter + 1] <- sum(C2newCof) # ClimVar2
      if(anova(Mod0, Mod)$RSS[1] > anova(Mod0, Mod)$RSS[2]){ # only save p value if model is an imporvement
        ps[counter+1] <- anova(Mod0, Mod)$'Pr(>F)'[2] 
      }else{ # if model is not an improvement over null, set p to 1
        ps[counter+1] <- 1}
      Mods_ls[[counter+1]] <- Mod # save model to list of models
      counter <- counter + 1 }
    ### Selecting best model -----
    AICs <- sapply(X = Mods_ls, FUN = AIC) # calculate AICs for each model
    # Best <- which(abs(AICs) == min(abs(AICs), na.rm = TRUE))[1] # best model, if same value present use first
    Best <- which(abs(coeffsC) == max(abs(coeffsC)))[1] # select fro strongest soil memory effect
    c_NDVI <- coeffst1[Best] # ndvi coefficient
    c_Clim <- coeffsC[Best] # climate coefficient
    c_Clim2 <- coeffsC2[Best] # climate 2 coefficient
    p_Mod <- ps[Best] # p-value is set to p-value of best model
    AICMod <- AICs[Best] # AIC is set to best AIC
    Bestlag <- Cumlags[Best] # this is the lag at which best model was observed
    ## EXPLAINED VARIANCE----
    colnames(ModData_df)[c(1:2)] <- c("NDVI_anom", "NDVI_Lag1")
    # Legendre & Legendre
    Explainedvar <- lm(data = ModData_df,
                       NDVI_anom ~ NDVI_Lag1 + ModData_df[,Best+5])
    Explainedvar <- summary(Explainedvar)[["r.squared"]]
    VarTotalNDVI <- lm(data = ModData_df,
                       NDVI_anom ~ NDVI_Lag1)
    VarTotalNDVI <- summary(VarTotalNDVI)[["r.squared"]]
    VarTotalQsoil <- lm(data = ModData_df,
                        NDVI_anom ~ ModData_df[,Best+5])
    VarTotalQsoil <- summary(VarTotalQsoil)[["r.squared"]]
    VarShared <- VarTotalQsoil + VarTotalNDVI - Explainedvar
    VarNDVI <- VarTotalNDVI - VarShared
    VarQsoil <- VarTotalQsoil - VarShared
    ## WRITING INFORMATION TO RASTERS----
    ModelEval_ras[pixel] <- as.numeric(c((Bestlag), AICMod, p_Mod, c_NDVI, c_Clim, c_Clim2, Explainedvar, 
                                         VarNDVI, VarShared, VarQsoil)) # saving model information to raster
    ## Updating progress bar----
    if(pixel == Data_Pos[1]){ # if we are currently on the first pixel
      T_End <- Sys.time() # note end time
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # calculate the time it took to establish and select models
      ## Put an estimator up on the console that tells the user when to expect the program to finish its current run
      print(paste("Calculating Vegetation Memory effects across ", Region, " should finish around: ", 
                  as.POSIXlt(T_Begin + Duration*length(Data_Pos), tz = Sys.timezone(location=FALSE)), sep=""))
      ## Update progress bar
      pb <- txtProgressBar(min = 0, max = length(Data_Pos), style = 3)
      pbi <- 0}
    pbi <- pbi + 1 ## Update progress bar
    setTxtProgressBar(pb, pbi)} # end of pixel loop
  ### Save raster ----- 
  writeRaster(ModelEval_ras, filename = paste(Dir.Memory,"/", Region, "_", ClimVar2, "-", ClimVar, 
                                              paste(Cumlags, collapse="_"),"_",FromY,"-", ToY, ".nc",sep=""),
              overwrite=TRUE, format="CDF")
  setwd(mainDir)}# end of VegMem function 

####--------------- CoeffScaling [ClimVar, ClimVar2, Region, Cumlags, FromY, ToY, UAbs] 
# (loading previously saved rasters and making visualisation of model coefficients better) ----
CoeffScaling <- function(ClimVar, ClimVar2, Region, Cumlags, FromY, ToY, UAbs){
  print("#################################################")
  print(paste("Producing composites of vegetation memory effects across ", unique(Region), sep=""))
  # PREPARATIONS ----
  Rasters <- ClimVar
  minmaxNDVIn <- rep(NA, length(ClimVar)*2)
  minmaxNDVIs <- rep(NA, length(ClimVar)*2)
  minmaxCVn <- rep(NA, length(ClimVar)*2)
  minmaxCVs <- rep(NA, length(ClimVar)*2)
  minmaxCV2n <- rep(NA, length(ClimVar)*2)
  minmaxCV2s <- rep(NA, length(ClimVar)*2)
  minmaxPos <- 1
  # LOADING DATA ----
  for(rasiter in 1:length(ClimVar)){ # cycle through specified vegetation memory rasters
    # load raster
    Alter_ras <- brick(paste(Dir.Memory,"/", Region[[rasiter]], "_", ClimVar2[[rasiter]], "-", ClimVar[[rasiter]], 
                             paste(Cumlags[[rasiter]], collapse="_"),"_",FromY,"-", ToY, ".nc",sep=""))
    Alter_ras <- Fun_NamesRas(raster = Alter_ras, ClimVar = ClimVar, ClimVar2 = ClimVar2, rasiter = rasiter)
    # PREPARING DATA ----
    P_ras <- Alter_ras$Model.p.value # extract p-value layer
    C_clim <- Alter_ras[[5]] # extract ClimVar coefficients
    C_climNon <- C_clim
    C_climNon[which(values(P_ras) < 0.05)] <- NA # set everything that significant to NA
    C_climSig <- C_clim
    C_climSig[which(values(P_ras) >= 0.05)] <- NA # set everything that's not significant to NA
    C_clim2 <- Alter_ras[[6]] # extract ClimVar2 coefficients
    C_climNon2 <- C_clim2
    C_climNon2[which(values(P_ras) < 0.05)] <- NA # set everything that significant to NA
    C_climSig2 <- C_clim2
    C_climSig2[which(values(P_ras) >= 0.05)] <- NA # set everything that's not significant to NA
    C_NDVI <- Alter_ras$Antecedent.NDVI..c_NDVI. # extract NDVI-1 coefficients
    C_NDVINon <- C_NDVI
    C_NDVINon[which(values(P_ras) < 0.05)] <- NA # set everything that significant to NA
    C_NDVISig <- C_NDVI
    C_NDVISig[which(values(P_ras) >= 0.05)] <- NA # set everything that's not significant to NA
    C_Lags <- Alter_ras[[1]] # extract Lags coefficients
    C_LagsNon <- C_Lags
    C_LagsNon[which(values(P_ras) < 0.05)] <- NA # set everything that significant to NA
    C_LagsSig <- C_Lags
    C_LagsSig[which(values(P_ras) >= 0.05)] <- NA # set everything that's not significant to NA
    # SAVING PARAMETERS ----
    Rasters[[rasiter]] <- list(C_NDVINon, C_NDVISig, C_climNon, C_climSig, 
                               C_climNon2, C_climSig2, C_LagsNon, C_LagsSig)
    ## Identify maximum values of each coefficient raster
    minmaxNDVIn[minmaxPos] <- max(values(C_NDVINon), na.rm = TRUE)
    minmaxNDVIs[minmaxPos] <-max(values(C_NDVISig), na.rm = TRUE)
    minmaxCVn[minmaxPos] <- max(values(C_climNon), na.rm = TRUE)
    minmaxCVs[minmaxPos] <- max(values(C_climSig), na.rm = TRUE)
    minmaxCV2n[minmaxPos] <- max(values(C_climNon2), na.rm = TRUE)
    minmaxCV2s[minmaxPos] <- max(values(C_climSig2), na.rm = TRUE)
    minmaxPos <- minmaxPos + 1 # +1 to counter in vector
    ## Identify minimum values of each coefficient raster
    minmaxNDVIn[minmaxPos] <- min(values(C_NDVINon), na.rm = TRUE)
    minmaxNDVIs[minmaxPos] <-min(values(C_NDVISig), na.rm = TRUE)
    minmaxCVn[minmaxPos] <- min(values(C_climNon), na.rm = TRUE)
    minmaxCVs[minmaxPos] <- min(values(C_climSig), na.rm = TRUE)
    minmaxCV2n[minmaxPos] <- min(values(C_climNon2), na.rm = TRUE)
    minmaxCV2s[minmaxPos] <- min(values(C_climSig2), na.rm = TRUE)
    minmaxPos <- minmaxPos + 1 # +1 to counter in vector
  } # end of loop cycling through specified vegetation memory rasters
  # MANN-WHITNEY U ----
  ## setting up directory
  Dir.Memory.Reg <- paste(Dir.Memory,"/",unique(Region),"-",FromY,"_",ToY, sep="")
  dir.create(Dir.Memory.Reg)
  ## cleaning directory o potential earlier runs
  if(paste("U-Variables_Abs",UAbs,".xlsx",sep="") %in% list.files(Dir.Memory.Reg)){
    file.remove(paste(Dir.Memory.Reg,"/U-Variables_Abs",UAbs,".xlsx",sep=""))}
  # Establish matrices and vectors for naming
  UModMat <- matrix(rep(NA, length(Rasters)^2), nrow=length(Rasters)) # for saving U outputs
  UVariables <- c("NDVI t-1", "Qsoil", "Tair", "Lags") # for naming purposes
  UMedians <- matrix(rep(NA, length(Rasters)^2), nrow=length(Rasters)) # for saving variable value medians
  dimnames(UMedians) <- list(c(1:4), UVariables) # set names
  # variable-wise comparison
  for(UVar in 1:length(Rasters)){ # loop over all variables
    for(UTest in 1:(length(Rasters)-1)){ # loop over the model layers
      UTest2 <- UTest + 1 # create seperate counter for variable with which to compare
      while(UTest2 <= length(Rasters)){ # cycle so long as second counter does not exceed range of specified models
        if(UAbs == TRUE){ # if absolute values should be used
          Test1 <- abs(values(Rasters[[UTest]][[UVar * 2]])) # data extraction
          Test2 <- abs(values(Rasters[[UTest2]][[UVar * 2]])) # data extraction
        }else{ # if absolute values are not desired
          Test1 <- values(Rasters[[UTest]][[UVar * 2]]) 
          Test2 <- values(Rasters[[UTest2]][[UVar * 2]])}
        test <- wilcox.test(Test1, Test2, paired = FALSE) # WHitney-U Test
        UModMat[UTest, UTest2] <- test$statistic # Extract test statistic
        UModMat[UTest2, UTest] <- test$p.value # Extract p-value
        Med1 <- median(Test1, na.rm=TRUE) # extract median
        Med2 <- median(Test2, na.rm=TRUE) # extract median
        UMedians[UTest,UVar] <- Med1 # write median of first object
        if(UTest2 == length(Rasters)){ # only write median of last variable
          UMedians[UTest2,UVar] <- Med2}# if statement
        UTest2 <- UTest2 + 1 } # while statement
    } # for UTest statement
    # saving output
    write.xlsx(UModMat, sheetName = UVariables[UVar], 
               file = paste(Dir.Memory.Reg,"/U-Variables_Abs",UAbs,".xlsx",sep=""), append = TRUE)
  } # for UVar statement
  write.xlsx(UMedians, sheetName = "Variable medians", # saving output
             file = paste(Dir.Memory.Reg,"/U-Variables_Abs",UAbs,".xlsx",sep=""), append = TRUE)
  UModelMat <- matrix(rep(NA,3^2),nrow=3) # matrix for model-internal comparisons
  dimnames(UModelMat) <- list(UVariables[1:3], UVariables[1:3]) # set names
  # model-wise comparison
  for(UTest in 1:(length(Rasters))){
    if(UAbs == TRUE){ # if absolute values should be used
      ND <- abs(values(Rasters[[UTest]][[2]]))
      QS <- abs(values(Rasters[[UTest]][[4]]))
      TA <- abs(values(Rasters[[UTest]][[6]]))
    }else{ # if absolute values are not desired
      ND <- values(Rasters[[UTest]][[2]])
      QS <- values(Rasters[[UTest]][[4]])
      TA <- values(Rasters[[UTest]][[6]])}
    UModelMat[2,1] <- wilcox.test(ND, QS, paired = FALSE)$p.value # WHitney-U Test
    UModelMat[3,2] <- wilcox.test(QS, TA, paired = FALSE)$p.value # WHitney-U Test
    UModelMat[3,1] <- wilcox.test(ND, TA, paired = FALSE)$p.value # WHitney-U Test
    UModelMat[1,2] <- wilcox.test(ND, QS, paired = FALSE)$statistic # WHitney-U Test
    UModelMat[2,3] <- wilcox.test(QS, TA, paired = FALSE)$statistic # WHitney-U Test
    UModelMat[1,3] <- wilcox.test(ND, TA, paired = FALSE)$statistic # WHitney-U Test
    write.xlsx(UModelMat, sheetName = paste("Model", UTest, sep=" "), # saving output
               file = paste(Dir.Memory.Reg,"/U-Variables_Abs",UAbs,".xlsx",sep=""), append = TRUE)
  } # for UTest statement
  # SAVING DATA FOR LATER PLOTTING ----
  for(iterplot in 1:length(ClimVar)){ # cycle through all specified vegetation memory raster for plotting
    ## lags -----
    Lag_ras <- brick(paste(Dir.Memory,"/", Region[[iterplot]], "_", ClimVar2[[iterplot]], "-", ClimVar[[iterplot]], 
                           paste(Cumlags[[iterplot]], collapse="_"),"_",FromY,"-", ToY, ".nc",sep=""))
    ## NDVI[t-1] -----
    Rasters[[iterplot]][[1]][1] <- max(minmaxNDVIn)
    Rasters[[iterplot]][[1]][2] <- min(minmaxNDVIn)
    Rasters[[iterplot]][[2]][3] <- max(minmaxNDVIs)
    Rasters[[iterplot]][[2]][4] <- min(minmaxNDVIs)
    ### ClimVar -----
    Rasters[[iterplot]][[3]][1] <- max(minmaxCVn)
    Rasters[[iterplot]][[3]][2] <- min(minmaxCVn)
    Rasters[[iterplot]][[4]][3] <- max(minmaxCVs)
    Rasters[[iterplot]][[4]][4] <- min(minmaxCVs)
    ### ClimVar2 -----
    Rasters[[iterplot]][[5]][1] <- max(minmaxCV2n)
    Rasters[[iterplot]][[5]][2] <- min(minmaxCV2n)
    Rasters[[iterplot]][[6]][3] <- max(minmaxCV2s)
    Rasters[[iterplot]][[6]][4] <- min(minmaxCV2s)
    ### sighnifciant coefficient rasters (already with max/min dots)
    Save_ras <- brick(Lag_ras[[1]], # Lags
                      Rasters[[iterplot]][[2]], # NDVI
                      Rasters[[iterplot]][[4]], # climvar
                      Rasters[[iterplot]][[6]]) # climvar2
    ### Saving significant effects -----
    writeRaster(Save_ras, filename = paste(Dir.Memory.Reg ,"/", ClimVar[[iterplot]], "_", ClimVar2[[iterplot]], 
                                           paste(Cumlags[[iterplot]], collapse="_"), "Plots.nc",sep=""),
                overwrite=TRUE, format="CDF")} # plotting loop
  setwd(mainDir)} # CoeffScaling end