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
  plot_ras <- crop(plot_ras, extent(Drylands))
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
  col.signega <- rev(got(n = 1000, alpha = 1, begin = 0, end = 1, direction = -1, option = "targaryen2"))
  col.sigposa <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "tyrell")
  col.signeg <- heat.colors(100)
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
  TriMap
  # COLOUR BARS, optional
  # ggframe <- data.frame(Values = na.omit(as.vector(Triplot)),
  #                       Idents = rep(SR_Titles, each = length(na.omit(as.vector(Triplot)))/3),
  #                       Order = rep(1:3, each = length(na.omit(as.vector(Triplot)))/3))
  # ggframeSmall <- na.omit(ggframe)
  # TriBars <- ggplotGrob(ggplot(ggframeSmall, aes(x = fct_reorder(Idents,Order), y = Values, fill = Idents)) + geom_boxplot() + theme_minimal(base_size = 15) + guides(fill=guide_legend(title="")) + theme(legend.position = c(.6, .9)) + ylab("Coefficients") +
  #                         theme(panel.background = element_rect(fill = "lightgrey")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  # )
  # GGCombined <- TriMap +
  #   annotation_custom(
  #     grob = TriBars,
  #     xmin = extent(plot_ras)[2]-5, # these coordinates need altering for different regions
  #     xmax = extent(plot_ras)[2],
  #     ymin = extent(plot_ras)[3],
  #     ymax = extent(plot_ras)[3]+5
  #   )
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
      plot(abs(plot_ras[[1]]), col=col.sigposa, legend=FALSE, axes=FALSE, add = TRUE)
      plot(abs(plot_ras[[1]]), legend.only=TRUE, col=col.sigposa, colNA = "black", 
           smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
    }else{
      if(Plot == 4 | Plot == 6){ # memory length
        plot(Countries, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
        plot(plot_ras[[Plot]], col = col.lags, axes=FALSE, add = TRUE, legend = FALSE)
        plot(abs(plot_ras[[Plot]]), legend.only=TRUE, col=col.lags, colNA = "black", 
             smallplot=smaplot, horizontal = TRUE, axis.args=list(cex.axis=1))
      }else{ # memory strength
        ## data
        Neg_ras <- plot_ras[[Plot]]
        Neg_ras[which(values(Neg_ras) >= 0)] <- NA
        Pos_ras <- plot_ras[[Plot]]
        Pos_ras[which(values(Pos_ras) < 0)] <- NA
        min <- min(values(plot_ras[[Plot]]), na.rm = TRUE)
        max <- max(values(plot_ras[[Plot]]), na.rm = TRUE)
        range <- abs(min) + max
        ## plotting
        plot(Countries, axes = FALSE, main = mainTit, col="#f2f2f2", bg="black", lwd=0.25)
        plot(Neg_ras, col=col.signeg, legend=FALSE, axes=FALSE, add = TRUE)
        plot(Pos_ras, col=col.sigpos, legend=FALSE, axes=FALSE, add = TRUE) 
        ## legend
        if(Scaled == FALSE){
          minsegment <- .805 * abs(min)/range
          smallplotxpos <- c(.13+minsegment,.935,.22,.25) # where to put colour scales
          smallplotxneg <- c(0.13,.13+minsegment,.22,.25) # where to put colour scales
        }else{
          minsegment <- .9 * abs(min)/range
          smallplotxpos <- c(.02+minsegment,.98,.17,.22) # where to put colour scales
          smallplotxneg <- c(0.02,.05+minsegment,.17,.22) # where to put colour scales
        }
        r <- raster(ncol=10, nrow=10) 
        values(r) <- seq(from = minValue(Neg_ras), to = 0, length=length(r))
        plot(r, legend.only=TRUE, col=col.signeg, colNA = "black", smallplot=smallplotxneg, horizontal = TRUE, axis.args=list(cex.axis=1))
        values(r) <- seq(from = 0, to = maxValue(Pos_ras), length=length(r))
        plot(r, legend.only=TRUE, col=col.sigpos, smallplot=smallplotxpos, horizontal = TRUE, axis.args=list(cex.axis=1))
      } # memory strength
    } # memory length
    if(Scaled == FALSE){
      dev.off()
    }
  } # plotting loop
  if(Scaled == TRUE){dev.off()}
  
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
  ## Plotting
  P_Diag <- ggplot(data = plot_df, aes(y = Data, x = Cell, fill = Variance)) + geom_bar(stat = "identity") + theme_bw(base_size= 25) + xlab("Raster Cells") + ylab("Variance") +
    scale_fill_manual(values=c(col.list[[1]], col.list[[2]], col.list[[3]], col.list[[4]],
                               col.list[[5]], col.list[[6]], col.list[[7]]))
  ggsave(plot = P_Diag, file=paste(Dir.Plots, "/", Region, "_VarParDiag.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)

  P_Box <- ggplot(data = plot_df, aes(y = Data, x = Variance, fill = Variance)) + geom_boxplot() + theme_bw(base_size= 25) + xlab("Raster Cells") + ylab("Variance") +
    scale_fill_manual(values=c(col.list[[1]], col.list[[2]], col.list[[3]], col.list[[4]],
                               col.list[[5]], col.list[[6]], col.list[[7]])) +
    theme(axis.text.x = element_text(size=7, angle=-30))
  ggsave(plot = P_Box, file=paste(Dir.Plots, "/", Region, "_VarParBox.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
  # P_VarsA <- P_Diag +
  #   annotation_custom(
  #     grob = ggplotGrob(P_Box + theme(legend.position = "none", 
  #                                     axis.title.x = element_blank(),
  #                                     axis.title.y = element_blank())),
  #     xmin = 0,
  #     xmax = ceiling(quantile(plot_df$Cell, .7)),
  #     ymin = .15,
  #     ymax = .5
  #   )
  # P_VarsB <- P_Box +
  #   annotation_custom(
  #     grob = ggplotGrob(P_Diag + theme(legend.position = "none", 
  #                                      axis.title.x = element_blank(),
  #                                      axis.title.y = element_blank())),
  #     xmin = 1.1,
  #     xmax = 6.9,
  #     ymin = .1,
  #     ymax = .25
  #   )
  
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
  if(Scaled == TRUE){
    Raster <- brick(paste(Dir.Memory, "/", Region, ".nc", sep=""))
    Raster[which(values(Raster[[15]]) > .05)] <- NA
    min <- min(mins)
    max <- max(maxs)
    col.mapview <- c(col.signeg[1:80], col.sigpos)
    col.breaks <- unique(c(seq(from = min, to = 0, length = 80), seq(from = 0, to = max, length = 100)))
    m0_c <- mapview(Countries, color = "black", alpha.regions = 0)
    m1_c <- mapview(layer.name = "Model AICs", 
                    abs(Raster[[1]]), legend = TRUE, col.regions = col.sigpos, maxpixels =  5755680, na.color = "#FFFFFF00")
    m2_c <- mapview(layer.name = "Intrinsic Memory (NDVI [t-1]",  at = col.breaks,
                    Raster[[2]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m3_c <- mapview(layer.name = "Qsoil1 Memory Effects",  at = col.breaks,
                    Raster[[3]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m4_c <- mapview(layer.name = "Qsoil1 Memory Length", at = 0:length(col.lags),
                    Raster[[4]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
    m5_c <- mapview(layer.name = "Tair Memory Effects",  at = col.breaks,
                    Raster[[5]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m6_c <- mapview(layer.name = "Tair Memory Length", at = 0:length(col.lags),
                    Raster[[6]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
  }else{
    col.mapview <- c(col.signeg[1:80], col.sigpos)
    m0_c <- mapview(Countries, color = "black", alpha.regions = 0)
    m1_c <- mapview(layer.name = "Model AICs", 
                    abs(Raster[[1]]), legend = TRUE, col.regions = col.sigpos, maxpixels =  5755680, na.color = "#FFFFFF00")
    m2_c <- mapview(layer.name = "Intrinsic Memory (NDVI [t-1]",  at = unique(c(seq(from = minValue(Raster[[2]]), to = 0, length = 80), seq(from = 0, to = maxValue(Raster[[2]]), length = 100))),
                    Raster[[2]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m3_c <- mapview(layer.name = "Qsoil1 Memory Effects",  at = unique(c(seq(from = minValue(Raster[[3]]), to = 0, length = 80), seq(from = 0, to = maxValue(Raster[[3]]), length = 100))),
                    Raster[[3]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m4_c <- mapview(layer.name = "Qsoil1 Memory Length", at = 0:length(col.lags),
                    Raster[[4]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
    m5_c <- mapview(layer.name = "Tair Memory Effects",  at = unique(c(seq(from = minValue(Raster[[5]]), to = 0, length = 80), seq(from = 0, to = maxValue(Raster[[5]]), length = 100))),
                    Raster[[5]], legend = TRUE, col.regions = col.mapview, maxpixels =  5755680, na.color = "#FFFFFF00")
    m6_c <- mapview(layer.name = "Tair Memory Length", at = 0:length(col.lags),
                    Raster[[6]], legend = TRUE, col.regions = col.lags, maxpixels =  5755680, na.color = "#FFFFFF00")
  }
  # Combine all maps into 1 map
  m_c <- m0_c + m1_c + m2_c + m3_c + m4_c+ m5_c+ m6_c
  # Plot the map
  # m_c
  mapshot(m_c, url = paste0(Dir.Plots, "/", Region,"_Effects_",Scaled,".html"))
  
  ## VARIANCE
  m0_c <- mapview(Countries, color = "black", alpha.regions = 0)
  m1_c <- mapview(layer.name = "Total Explained Variance", at = seq(from = 0, to = maxValue(Alter_ras[[1]]), length = 100), 
                  Alter_ras[[1]], legend = TRUE, col.regions = col.varpar1, maxpixels =  5755680, na.color = "#FFFFFF00")
  m2_c <- mapview(layer.name = "NDVI[t-1]", at = seq(from = 0, to = maxValue(Alter_ras[[2]]), length = 100),
                  Alter_ras[[2]], legend = TRUE, col.regions = col.varpar1, maxpixels =  5755680, na.color = "#FFFFFF00")
  m3_c <- mapview(layer.name = "Qsoil1",  at = seq(from = 0, to = maxValue(Alter_ras[[3]]), length = 100),
                  Alter_ras[[3]], legend = TRUE, col.regions = col.varpar2, maxpixels =  5755680, na.color = "#FFFFFF00")
  m4_c <- mapview(layer.name = "Tair",  at = seq(from = 0, to = maxValue(Alter_ras[[4]]), length = 100),
                  Alter_ras[[4]], legend = TRUE, col.regions = col.varpar3, maxpixels =  5755680, na.color = "#FFFFFF00")
  m5_c <- mapview(layer.name = "NDVI[t-1] + Qsoil1",  at = seq(from = 0, to = maxValue(Alter_ras[[5]]), length = 100),
                  Alter_ras[[5]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m6_c <- mapview(layer.name = "NDVI[t-1] + Tair",  at = seq(from = 0, to = maxValue(Alter_ras[[6]]), length = 100),
                  Alter_ras[[6]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m7_c <- mapview(layer.name = "Qsoil1 + Tair",  at = seq(from = 0, to = maxValue(Alter_ras[[7]]), length = 100),
                  Alter_ras[[7]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  m8_c <- mapview(layer.name = "Shared by all",  at = seq(from = 0, to = maxValue(Alter_ras[[8]]), length = 100),
                  Alter_ras[[8]], legend = TRUE, col.regions = col.varpar4, maxpixels =  5755680, na.color = "#FFFFFF00")
  # Combine all maps into 1 map
  m_c <- m0_c+ m1_c + m2_c + m3_c + m4_c + m5_c + m6_c + m7_c + m8_c
  # Plot the map
  # m_c
  mapshot(m_c, url = paste0(Dir.Plots, "/", Region,"_Variances.html"))
} # end of Fun_Plot

################################################################################
####--------------- Fun_Region [] # assess correllations of different vegmem characteristics in different regions)
Fun_PlotReg <- function(){
  ### DATA ----
  Raster <- brick(paste(Dir.Memory, "/GlobalDrylands.nc", sep=""))
  ### Ecoregions ---
  if(!file.exists(file.path(Dir.Mask, "WWF_ecoregions"))){
    # downloading Terrestrial Ecoregion Shapefile as zip
    download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                  destfile = file.path(Dir.Mask, "wwf_ecoregions.zip")
    )
    # unpacking the zip
    unzip(file.path(Dir.Mask, "wwf_ecoregions.zip"),
          exdir = file.path(Dir.Mask, "WWF_ecoregions")
    )
  }
  wwf <- readOGR(file.path(Dir.Mask, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"), verbose = FALSE) # loading shapefile for biomes
  # wwf <- rgeos::gBuffer(wwf, byid=TRUE, width=0)
  # wwf <- crop(wwf, extent(Raster)) # crop to extent of vegmem
  
  ### Realms ---
  PlotColsReg <- rainbow(length(unique(wwf@data[["REALM"]])))
  wwf@data[,"REALM"] <- as.factor(wwf@data[,"REALM"])
  parcel_ras <- rasterize(x = wwf, y = Raster, field = "REALM", fun = 'first') # rasterize
  levelplot(parcel_ras, at = c(sort(unique(values(parcel_ras))), max(unique(values(parcel_ras)), na.rm = TRUE)+1) -.5, 
            col.regions = PlotColsReg)
  
  ### Biomes ---
  parcel_ras2 <- rasterize(x = wwf, y = Raster, field = "BIOME", fun = 'first') # rasterize
  values(parcel_ras2)[which(values(parcel_ras2) == 98 | values(parcel_ras2) == 99)] <- NA
  PlotColsReg2 <- rainbow(length(unique(values(parcel_ras2))))
  levelplot(parcel_ras2, at = c(sort(unique(values(parcel_ras2))), max(unique(values(parcel_ras2)), na.rm = TRUE)+1) -.5, 
            col.regions = PlotColsReg2)
  
  ### Raster stacking ----
  Raster <- stack(Raster, parcel_ras, parcel_ras2) # combine raster
  Raster_df <- na.omit(as.data.frame(Raster)) # convert to data frame and omit NAs
  colnames(Raster_df) <- c("Model AICs", "NDVI [t-1]", "Qsoil1", "Lag Qsoil1", "Tair", "Lag Tair",
                           "Total Explained Variance", "V_NDVI", "V_C1", "V_C2", "V_C1.2", 
                           "V_C1.3", "V_C2.3", "V_Shared", "Model p-value", "REALM", "BIOME") # set column names
  ### Recoding Realms ---
  recode <- paste0(sort(unique(Raster_df$REALM)), 
                   paste0(" = '", levels(wwf@data[["REALM"]])[sort(unique(Raster_df$REALM))], "'"), 
                   collapse = ";")
  Raster_df$REALM <- recode(Raster_df$REALM,recode)
  
  ### Establishing Biomes in Realms ---
  Raster_df$RELBIO <- paste(Raster_df$REALM, Raster_df$BIOME, sep ="_") # combine biomes and realms
  Raster_df$RELBIO[which(Raster_df$RELBIO %in% names(which(table(Raster_df$RELBIO) < 50)))] <- NA # omit biomes in realms with less than 50 samples
  Raster_df <- na.omit(Raster_df)
  
  ### RECODING ABBREVIATIONS ---
  Abbr_Realms <- levels(wwf@data[["REALM"]])
  Full_Realms <- c("Australasia", "Antarctic", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")
 for(Iter_Realms in 1:length(Abbr_Realms)){
   Raster_df$REALM[which(Raster_df$REALM == Abbr_Realms[Iter_Realms])] <- Full_Realms[Iter_Realms]
 }
  
  Abbr_Biomes <- 1:18
  Full_Biomes <- c("Tropical & Subtropical Moist Broadleaf Forests",
                   "Tropical & Subtropical Dry Broadleaf Forests",
                   "Tropical & Subtropical Coniferous Forests",
                   "Temperate Broadleaf & Mixed Forests",
                   "Temperate Conifer Forests",
                   "Boreal Forests/Taiga",
                   "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                   "Temperate Grasslands, Savannas & Shrublands",
                   "Flooded Grasslands & Savannas",
                   "Montane Grasslands & Shrublands",
                   "Tundra",
                   "Mediterranean Forests, Woodlands & Scrub",
                   "Deserts & Xeric Shrublands",
                   "Mangroves")
  for(Iter_Biomes in 1:length(Abbr_Biomes)){
    Raster_df$BIOME[which(Raster_df$BIOME == Abbr_Biomes[Iter_Biomes])] <- Full_Biomes[Iter_Biomes]
  }
  
  ### TOTAL EXPLAINED VARIANCE ~ INTRINSIC ----
  X = "NDVI [t-1]"
  Y = "Total Explained Variance"
  Z = "Tair"
  ## limit data frame to desired variables
  Raster_df2 <- data.frame(X = rep(Raster_df[, which(colnames(Raster_df) == Y)], 2),
                          Y = c(Raster_df[, which(colnames(Raster_df) == X)], Raster_df[, which(colnames(Raster_df) == Z)]),
                          Effect = rep(c(X, Z), each = length(Raster_df[, which(colnames(Raster_df) == Y)])),
                          REALM = rep(Raster_df$REALM, 2),
                          BIOME = rep(Raster_df$BIOME, 2),
                          RELBIO = rep(Raster_df$RELBIO, 2))
  
  ### PLOTTING
  p_REALMS <- ggplot(Raster_df2, aes(x = Y, y = X, linetype = factor(REALM), color = as.factor(Effect))) + 
    stat_smooth(method="lm", se = TRUE) + 
    labs(x = "Effect Size", y = Y) +
    labs(linetype='Realm', color = "Memory Effect") + 
    ylim(c(0, 1)) + guides(linetype = FALSE) +
    theme_bw() + scale_colour_manual(values = c("darkgreen", "red"))
  
  
  ## limit data frame to desired variables
  Raster_df2 <- data.frame(X = Raster_df[, which(colnames(Raster_df) == X)],
                           Y = Raster_df[, which(colnames(Raster_df) == Y)],
                           REALM = Raster_df$REALM,
                           BIOME = Raster_df$BIOME,
                           RELBIO = Raster_df$RELBIO)
  ### REALMS ---
  NDVI_null <- glm(Y ~ 1, data = Raster_df2)
  NDVI_glm <- glm(Y ~ X*REALM, data = Raster_df2)
  anova(NDVI_null, NDVI_glm)
  ### BIOMES in REALMS ---
  NDVI_null <- glm(Y ~ 1, data = Raster_df2)
  NDVI_glm <- glm(Y ~ X*RELBIO, data = Raster_df2)
  anova(NDVI_null, NDVI_glm)
  ### PLOTTING
  p_RELBIO <- ggplot(Raster_df2, aes(x = X, y = Y, color = factor(BIOME), linetype = factor(REALM))) + 
    stat_smooth(method="lm", se = TRUE) + 
    labs(x = X, y = Y) +
    labs(color='Biome', linetype = "Realm") + 
    ylim(c(0, 1)) +
    theme_bw() + scale_color_viridis(discrete = TRUE, option = "D") +
    guides(color=guide_legend(nrow = 5), linetype=guide_legend(nrow = 6))
  
  prow <- plot_grid(p_REALMS + theme(legend.position='none'), 
                    p_RELBIO + theme(legend.position='none'), 
                    align = 'vh',
                    labels = c("A", "B"),
                    hjust = -1,
                    nrow = 1)
  legend_a <- get_legend(p_REALMS + theme(legend.position="bottom"))
  legend_b <- get_legend(p_RELBIO + theme(legend.position="bottom"))
  p <- plot_grid(prow, legend_a, legend_b, ncol = 1, rel_heights = c(1, .07, .35))
  ggsave(p, file=paste(Dir.Plots, "/RegionalDifferencesTEST.jpeg", sep = ""), width = 32, height = 21, units = "cm", quality = 100)
  
  ### INTRINSIC ~ LAG SOIL ----
  X = "Lag Qsoil1"
  Y = "NDVI [t-1]"
  ## limit data frame to desired variables
  Raster_df2 <- data.frame(X = Raster_df[, which(colnames(Raster_df) == X)],
                           Y = Raster_df[, which(colnames(Raster_df) == Y)],
                           REALM = Raster_df$REALM,
                           RELBIO = Raster_df$RELBIO,
                           BIOME = Raster_df$BIOME)
  
  ### REALMS ---
  NDVI_null <- glm(Y ~ 1, data = Raster_df2)
  NDVI_glm <- glm(Y ~ X*REALM, data = Raster_df2)
  anova(NDVI_null, NDVI_glm)
  ### PLOTTING
  p1_REALMS <- ggplot(Raster_df2, aes(x = X, y = Y, linetype = factor(REALM), linewidth = 2)) + 
    stat_smooth(method="lm", se = TRUE, size = 1.2) + 
    labs(x = X, y = Y) +
    labs(linetype='Realm') +
    theme_bw()
  
  ### BIOMES in REALMS ---
  Raster_df3 <- Raster_df2[which(Raster_df2$REALM == "Australasia"), ] # limiting to just one realm
 
  NDVI_null <- glm(Y ~ 1, data = Raster_df3)
  NDVI_glm <- glm(Y ~ X*BIOME, data = Raster_df3)
  anova(NDVI_null, NDVI_glm)
  ### PLOTTING
  p1_RELBIO <- ggplot(Raster_df3, aes(x = X, y = Y, color = factor(BIOME))) +
    stat_smooth(method="lm", se = TRUE) + 
    labs(x = X, y = Y) +
    labs(color='Biomes in Australasia') + 
    theme_bw()  + scale_color_viridis(discrete = TRUE, option = "D")
  
  g <- plot_grid(p1_REALMS, 
                    p1_RELBIO, 
                    align = 'vh',
                    labels = c("A", "B"),
                    hjust = -1,
                    nrow = 2)
  ggsave(g, file=paste(Dir.Plots, "/RegionalDifferences2.jpeg", sep = ""), width = 32, height = 18, units = "cm", quality = 100)
} # end of function