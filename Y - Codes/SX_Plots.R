################################################################################
####--------------- Fun_Plot [Region, Scaled] # all things vegetation memory (maps, varpar, system dynamics)
Fun_Plot <- function(Region, Scaled = FALSE){
  ####--------------- FILE SELECTION ----------------
  Raster <- brick(paste(Dir.Memory, "/", Region, ".nc", sep=""))
  Raster[which(values(Raster[[15]]) > .05)] <- NA
  Back <- raster(paste(Dir.KrigCov, "Co-variates_NativeResolution.nc", sep="/"),
                 varname = "Elevation")
  Countries <- readOGR(Dir.Mask, "ne_50m_admin_0_countries", verbose = FALSE)
  setwd(Dir.Mask)
  Drylands <- shapefile("dryland_2")
  Raster <- crop(Raster, extent(Drylands))
  Back <- crop(Back, extent(Drylands))
  Countries <- crop(Countries, extent(Drylands))
  ####--------------- MISC ---------------- 
  SR_Titles <- list("NDVI[t-1]", "Air Temperature", "Soil Moisture (0-7cm)")
  height <- 11
  width <- 16
  ####--------------- SCALING -----------------
  if(Scaled == TRUE){
    # loop over all different model coefficients
    counter = 1
    maxs <- NA
    mins <- NA
    for(rasterscale in c(2,3,5)){
      maxs[counter] <- maxValue(Raster[[rasterscale]])
      mins[counter] <- minValue(Raster[[rasterscale]])
      counter <- counter + 1
    }
    # fill first and second NA with max and min
    for(rasterscale in c(2,3,5)){
      values(Raster[[rasterscale]])[which(!is.na(values(Raster[[rasterscale]])))[97:100]] <- max(maxs)
      values(Raster[[rasterscale]])[which(!is.na(values(Raster[[rasterscale]])))[100:105]] <- min(mins)
    }
  }
  
  ####--------------- PLOTTING ----------------
  ##------- MEMORY COMPONENTS -------
  ## Plotting Setup
  col.signega <- rev(got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "targaryen2"))
  col.sigposa <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "tyrell")
  col.signeg <- heat.colors(100)[1:70]
  col.sigpos <- rev(viridis(100))
  col.nonsig <- colorRampPalette(c("grey"))(1)
  col.lags <- got(n = 12, alpha = 1, begin = 0, end = 1, direction = 1, option = "daenerys")
  ##------- TRICOLOUR -------
  Tricols <- list()
  plot_ras <- Raster[[c(2,3,5)]]
  Triplot <- abs(plot_ras)
  TriMap <- ggR(Back) + 
    ggRGB(img = Triplot, r = 3, g = 1, b = 2, scale = max(maxValue(Triplot)), 
          stretch = 'none', ggLayer = TRUE) +
    ggtitle(label = "Relative Importance of Vegetation Responses") + 
    ylab("Latitude [°]") + xlab("Longitude [°]") + 
    theme_bw(base_size = 15)
  ggsave(file=paste(Dir.Plots, "/", Region, "_RelImportance.jpeg", sep = ""), width = width, height = height, units = "cm", quality = 100)
  
  ##------- MODEL COEFFICIENTS -------
  plot_ras <- Raster[[c(1:6)]]
  ### COEFFICIENTS
  IndTitles <- c("Model AICs",
                 paste(SR_Titles[1], "(Intrinsic Memory)"),
                 paste(SR_Titles[3], "(Inverse Resistance)"),
                 paste(SR_Titles[3], "(Memory Length)"),
                 paste(SR_Titles[2], "(Inverse Resistance)"),
                 paste(SR_Titles[2], "(Memory Length)"))
  SaveTitels <- c("AIC", "S_NDVI", "S_Qsoil1", "L_Qsoil1", "S_Tair", "L_Tair")
  
  if(Scaled == TRUE){
    jpeg(file=paste(Dir.Plots, "/", Region, "_SCALEDOverview.jpeg", sep = ""), width = width*1.5, height = height*2, units = "cm", quality = 100, res = 1000)
    par(mfrow = c(3,2), mai = c(.7, 0, 0.2, 0))
    smaplot <- c(.02, .98, .17, .22)
  }
  for(Plot in c(1,2,4,3,6,5)){
    mainTit <- IndTitles[Plot]
    if(Scaled == FALSE){
      jpeg(file=paste(Dir.Plots, "/", Region, "_", SaveTitels[Plot],".jpeg", sep = ""), 
           width = width, height = height, units = "cm", quality = 100, res = 1000)
      smaplot <- c(.13, .935, .22, .25)
    }
    if(Plot == 1){
      plot(Countries, axes = FALSE, main = "Model AICs", col="#f2f2f2", bg="black", lwd=0.25)
      axis(1, at=seq(from = -180, to = 180, by = 20), 
           labels = paste(seq(from = -180, to = 180, by = 20), "°"), 
           tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
      axis(2, at=seq(from = -50, to = 60, by = 10), 
           labels = paste(seq(from = -50, to = 60, by = 10), "°"), 
           tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
      plot(abs(plot_ras[[1]]), col=col.sigposa, legend=FALSE, axes=FALSE, add = TRUE)
      plot(abs(plot_ras[[1]]), legend.only=TRUE, col=col.sigposa, colNA = "black", 
           smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
    }else{
      if(Plot == 4 | Plot == 6){ # memory length
        plot(Countries, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
        axis(1, at=seq(from = -180, to = 180, by = 20), 
             labels = paste(seq(from = -180, to = 180, by = 20), "°"), 
             tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
        axis(2, at=seq(from = -50, to = 60, by = 10), 
             labels = paste(seq(from = -50, to = 60, by = 10), "°"), 
             tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
        plot(plot_ras[[Plot]], col = col.lags, axes=FALSE, add = TRUE, legend = FALSE)
        plot(abs(plot_ras[[Plot]]), legend.only=TRUE, col=col.lags, colNA = "black", 
             smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
      }else{ # memory strength
        plot(Countries, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
        axis(1, at=seq(from = -180, to = 180, by = 20), 
             labels = paste(seq(from = -180, to = 180, by = 20), "°"), 
             tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
        axis(2, at=seq(from = -50, to = 60, by = 10), 
             labels = paste(seq(from = -50, to = 60, by = 10), "°"), 
             tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
        plot(plot_ras[[Plot]], col=col.sigpos, legend=FALSE, axes=FALSE, add = TRUE)
        plot(plot_ras[[Plot]], legend.only=TRUE, col=col.sigpos, colNA = "black", 
             smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
      } # memory strength
    } # memory length
    if(Scaled == FALSE){
      dev.off()
    }
  } # plotting loop
  if(Scaled == TRUE){dev.off()}
  
  if(Scaled == TRUE){
    AUS_shp <- Countries[which(Countries$NAME_SORT == "Australia"),]
    AUS_ras <- crop(plot_ras, extent(c(extent(AUS_shp)[1], extent(AUS_shp)[2], -40, -12)))
    Countries_AU <- crop(Countries, extent(c(extent(AUS_shp)[1], extent(AUS_shp)[2], -40, -12)))
    jpeg(file=paste(Dir.Plots, "/", Region, "_SCALEDAUS.jpeg", sep = ""), width = width*1.5, height = height, units = "cm", quality = 100, res = 1000)
    par(mfrow = c(1,3), mai = c(.7, 0, 0.2, 0))
    smaplot <- c(.02, .98, .1, .15)
    for(Plot in c(2,3,4)){
      mainTit <- IndTitles[Plot]
        if(Plot == 4){ # memory length
          plot(Countries_AU, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
          axis(1, at=seq(from = 80, to = 180, by = 10), 
               labels = paste(seq(from = 80, to = 180, by = 10), "°"), 
               tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
          axis(2, at=seq(from = -35, to = -15, by = 5), 
               labels = paste(seq(from = -35, to = -15, by = 5), "°"), 
               tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
          plot(AUS_ras[[Plot]], col = col.lags, axes=FALSE, add = TRUE, legend = FALSE)
          plot(abs(AUS_ras[[Plot]]), legend.only=TRUE, col=col.lags, colNA = "black", 
               smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
        }else{ # memory strength
          plot(Countries_AU, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
          axis(1, at=seq(from = 80, to = 180, by = 10), 
               labels = paste(seq(from = 80, to = 180, by = 10), "°"), 
               tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
          axis(2, at=seq(from = -35, to = -15, by = 5), 
               labels = paste(seq(from = -35, to = -15, by = 5), "°"), 
               tck=.02, col = "white", col.axis = "white", mgp = c(0,-1.4,0), cex.axis = .7)
          plot(AUS_ras[[Plot]], col=col.sigpos, legend=FALSE, axes=FALSE, add = TRUE)
          plot(AUS_ras[[Plot]], legend.only=TRUE, col=col.sigpos, colNA = "black", 
               smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
        } # memory strength
      } # memory length
    dev.off()
    } # plotting loop
  
  # ##------- VARIANCE PARTITIONING -------
  ## Plotting Setup
  col.list <- as.list(c(got(n = 3, alpha = 1, begin = 0.2, end = 1, direction = 1, option = "tyrell"),
                        got(n = 2, alpha = 1, begin = 0.2, end = 1, direction = 1, option = "white_walkers"),
                        got(n = 2, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "margaery")[1],
                        got(n = 2, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen2")[1]))
  ## Data
  Alter_ras <- Raster[[7:14]]
  Alter_ras[2] <- 0
  values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
  values(Alter_ras)[which(values(Alter_ras) > 1)] <- 1
  Order_ras <- sum(Alter_ras[[2:8]], na.rm = FALSE)
  cells <- order(values(Order_ras), na.last = TRUE) # ordering and putting NAs last
  cells <- cells[1:(which(is.na(values(Order_ras)[cells]))[1]-1)] # only keeping cells before first NA
  plot_df <- data.frame(Data =NA, Cell =NA, Variance =NA)
  Idents <- c("Total", "NDVI[t-1]", "Qsoil1", "Tair", 
              "NDVI[t-1] + Qsoil1", "NDVI[t-1] + Tair", "Qsoil1 + Tair", "Shared by all")
  for(i in 2:8){
    values(Alter_ras[[i]])[cells][which(values(Alter_ras[[i]])[cells] > quantile(values(Alter_ras[[i]])[cells], .99))] <- quantile(values(Alter_ras[[i]])[cells], .95)
    plot_df1 <- data.frame(Data = values(Alter_ras[[i]])[cells],
                           Cell = c(1:length(cells)),
                           Variance = rep(Idents[i], length(cells)))
    plot_df <- rbind(plot_df, plot_df1)
  }
  plot_df <- na.omit(plot_df)
  P_Box <- ggplot(data = plot_df, aes(y = Data, x = Variance, fill = Variance)) + geom_boxplot() + theme_bw(base_size= 25) + xlab("Variable Components") + ylab("Proportion Variance Explained") +
    scale_fill_manual(values=c(col.list[[1]], col.list[[2]], col.list[[3]], col.list[[4]],
                               col.list[[5]], col.list[[6]], col.list[[7]])) +
    theme(axis.text.x = element_text(size=20, angle=-30)) + theme(legend.position = "none")
  ggsave(plot = P_Box, file=paste(Dir.Plots, "/", Region, "_VarParBox.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
  
  # Maps
  col.varpar1 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "wildfire")
  col.varpar2 <- got(n = 100, alpha = 1, begin = .4, end = 1, direction = 1, option = "jon_snow")
  col.varpar3 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen")
  col.varpar4 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "margaery")
  col.list <- list(col.varpar1, col.varpar1, col.varpar2, col.varpar3,
                   col.varpar4, col.varpar4, col.varpar4, col.varpar4)
  Titles <- c("Total Explained Variance", "NDVI[t-1]", "Qsoil1", "Tair", "NDVI[t-1] + Qsoil1",
              "NDVI[t-1] + Tair", "Qsoil1 + Tair", "Shared by all")
  ## singular
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParMap1.jpeg", sep = ""), width = 32, height = 17, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(2,2))
  Alter_ras <- Raster[[7:14]]
  Alter_ras[2] <- 0
  values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
  for(i in 1:4){
    col.varpar <- col.list[[i]]
    plot(Alter_ras[[i]], col=col.varpar, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[i])
    plot(Alter_ras[[i]], legend.only=TRUE, col=col.varpar, smallplot=c(.1, .9, .2, .23), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  }
  dev.off()
  ## shared
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParMap2.jpeg", sep = ""), width = 32, height = 17, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(2,2))
  Alter_ras <- Raster[[7:14]]
  Alter_ras[2] <- 0
  values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
  for(i in 5:8){
    col.varpar <- col.list[[i]]
    plot(Alter_ras[[i]], col=col.varpar, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[i])
    plot(Alter_ras[[i]], legend.only=TRUE, col=col.varpar, smallplot=c(.1, .9, .2, .23), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  }
  dev.off()
  
  # ##------- MAPVIEW -------
  ## EFFECTS
  Raster <- brick(paste(Dir.Memory, "/", Region, ".nc", sep=""))
    Raster[which(values(Raster[[15]]) > .05)] <- NA
    col.mapview <- col.sigpos
    m0_c <- mapview(Countries, color = "black", alpha.regions = 0)
    # m1_c <- mapview(layer.name = "Model AICs", 
    #                 abs(Raster[[1]]), legend = TRUE, col.regions = col.sigpos, maxpixels =  5755680, na.color = "#FFFFFF00")
    m2_c <- mapview(layer.name = "Intrinsic Memory (NDVI [t-1]",
                    Raster[[2]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m3_c <- mapview(layer.name = "Qsoil1 Memory Effects",
                    Raster[[3]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m4_c <- mapview(layer.name = "Qsoil1 Memory Length", at = 0:length(col.lags),
                    Raster[[4]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
    m5_c <- mapview(layer.name = "Tair Memory Effects",
                    Raster[[5]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m6_c <- mapview(layer.name = "Tair Memory Length", at = 0:length(col.lags),
                    Raster[[6]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
  # Combine all maps into 1 map
  m_c <- m0_c + m2_c + m3_c + m4_c+ m5_c+ m6_c
  # Plot the map
  # m_c
  mapshot(m_c, url = paste0(Dir.Plots, "/", Region,"_Effects.html"))
  
  ## VARIANCE
  m0_c <- mapview(Countries, color = "black", alpha.regions = 0)
  m1_c <- mapview(layer.name = "Total Explained Variance",
                  Alter_ras[[1]], legend = TRUE, col.regions = col.varpar1, maxpixels =  5755680, na.color = "#FFFFFF00")
  m2_c <- mapview(layer.name = "NDVI[t-1]",
                  Alter_ras[[2]], legend = TRUE, col.regions = col.varpar1, maxpixels =  5755680, na.color = "#FFFFFF00")
  m3_c <- mapview(layer.name = "Qsoil1",
                  Alter_ras[[3]], legend = TRUE, col.regions = col.varpar2, maxpixels =  5755680, na.color = "#FFFFFF00")
  m4_c <- mapview(layer.name = "Tair",
                  Alter_ras[[4]], legend = TRUE, col.regions = col.varpar3, maxpixels =  5755680, na.color = "#FFFFFF00")
  m5_c <- mapview(layer.name = "NDVI[t-1] + Qsoil1",
                  Alter_ras[[5]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m6_c <- mapview(layer.name = "NDVI[t-1] + Tair",
                  Alter_ras[[6]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m7_c <- mapview(layer.name = "Qsoil1 + Tair",
                  Alter_ras[[7]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m8_c <- mapview(layer.name = "Shared by all",
                  Alter_ras[[8]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  # Combine all maps into 1 map
  m_c <- m0_c+ m1_c + m2_c + m3_c + m4_c + m5_c + m6_c + m7_c + m8_c
  # Plot the map
  # m_c
  mapshot(m_c, url = paste0(Dir.Plots, "/", Region,"_Variances.html"))
} # end of Fun_Plot