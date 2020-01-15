################################################################################
####--------------- Fun_Plot [Region, SoilLayer] # all things vegetation memory (maps, varpar, system dynamics)
Fun_Plot <- function(Region, SoilLayer = 1, Scaled = FALSE){
  Region <- unlist(Region)
  ####--------------- FILE SELECTION ----------------
  Rasters <- list()
  for(i in 1:length(Region)){
    Files <- list.files(Dir.Memory)[grep(pattern = Region[i], list.files(Dir.Memory))]
    Files <- Files[grep(pattern = paste("Qsoil", SoilLayer, sep=""), Files)]
    File <- Files[grep(pattern = ".nc", Files)]
    Memory <- brick(paste(Dir.Memory,File, sep="/"))
    Rasters[[i]] <- Memory
  }
  ####--------------- MISC ----------------
  SR_Titles <- list("NDVI[t-1]", "Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", "Soil Moisture (28-100cm)", "Soil Moisture (100-255cm)")
  Titles <- c(paste(SR_Titles[1], "(Intrinsic Memory)"), paste(SR_Titles[2+SoilLayer], "(Inverse Resistance)"), paste(SR_Titles[2], "(Inverse Resistance)"))
  TitlesShort <- c(SR_Titles[1][[1]], SR_Titles[2+SoilLayer][[1]], SR_Titles[2][[1]])
  ####--------------- SCALING -----------------
  if(Scaled == TRUE){
    # loop over all different model coefficients
    for(rasterscale in c(1,4:6)){
      maxs <- NA
      mins <- NA
      # identify maxs and mins
      for(scale in 1:length(Region)){
        maxs[scale] <- maxValue(Rasters[[scale]][[rasterscale]])
        mins[scale] <- minValue(Rasters[[scale]][[rasterscale]])
      }
      # fill first and second NA with max and min
      for(scale in 1:length(Region)){
        values(Rasters[[scale]][[rasterscale]])[which(is.na(values(Rasters[[scale]][[rasterscale]])))[1]] <- max(maxs)
        values(Rasters[[scale]][[rasterscale]])[which(is.na(values(Rasters[[scale]][[rasterscale]])))[2]] <- min(mins)
      }
    }
  }
  ####--------------- PLOTTING ----------------
  ##------- MEMORY COMPONENTS -------
  ## Plotting Setup
  # col.signeg <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "targaryen2")
  # col.sigpos <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "tyrell")
  col.signeg <- heat.colors(100)
  col.sigpos <- rev(viridis(100))
  col.nonsig <- colorRampPalette(c("grey"))(1)
  col.lags <- got(n = 12, alpha = 1, begin = 0, end = 1, direction = 1, option = "daenerys")
  ##------- TRICOLOUR -------
  Tricols <- list()
  for(plot in 1:length(Region)){
    plot_ras <- Rasters[[plot]][[c(1,4:6)]]
    Triplot <- abs(plot_ras[[2:4]])
    TriMap <- ggRGB(img = Triplot, r = 3, g = 1, b = 2, scale = max(maxValue(Triplot)), stretch = 'none') + ggtitle(label = "Relative Importance of Vegetation Responses") + ylab("Latitude [°]") + xlab("Longitude [°]") + theme_bw(base_size = 15)
    Tricols[[plot]] <- TriMap
  
    ## COLOUR BARS 
    # ggframe <- data.frame(Values = as.vector(Triplot),
    #                       Idents = rep(TitlesShort, each = dim(values(Triplot))[1]),
    #                       Order = rep(1:3, each = dim(values(Triplot))[1]))
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
    ggsave(file=paste(Dir.Plots, "/", Region[plot], "_RelImportance", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
  }
  
  ##------- MODEL COEFFICIENTS -------
  ## setting jpeg size
  if(length(Region) == 1){
    height <- 22
    width <- 32
  }else{
    height <- 11*length(Region)
    width <- 16*2
  }
  
  ## plotting to jpeg
  jpeg(file=paste(Dir.Plots, "/", paste(Region,collapse=""), "_Model", SoilLayer, "_", Scaled, ".jpeg", sep = ""), width = width, height = height, units = "cm", quality = 100, res = 1000)
  if(length(Region)>1){
    if(Scaled == TRUE){
      par(mfrow=c(length(Region), 4),
          oma = c(3,0,1,0) + 0.1,
          mar = c(3,.25,3,0) + 0.1)
    }else{
      par(mfrow=c(length(Region), 4),
          oma = c(3,0,1,0) + 0.1,
          mar = c(3,.25,3,0) + 0.1) # margins: bottom, left, top and right 
    }
  }else{
    par(mfrow=c(2,2))
  }
  for(plot in 1:length(Region)){
    plot_ras <- Rasters[[plot]][[c(1,4:6)]]
    if(Scaled == FALSE | Scaled == TRUE & plot == 1){
      mainTit <- paste(SR_Titles[2+SoilLayer][[1]], "(Memory Length)")
    }else{
      mainTit <- ""
    }
    plot(plot_ras[[1]], col = col.lags, colNA = "black", legend=FALSE, axes=FALSE, main = mainTit)
    if(length(Region) > 1 & Scaled == FALSE | Scaled == TRUE & plot == length(Region)){ # plot legend only for last raster stack if we scale
      plot(plot_ras[[1]], legend.only=TRUE, col=col.lags, smallplot=c(0.017, .93, .08, .105), horizontal = TRUE, axis.args=list(cex.axis=1.5))
    }else{
      if(Scaled == FALSE){
        plot(plot_ras[[1]], legend.only=TRUE, col=col.lags, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=1.5))  
      }
    }
    for(Plot in 2:4){
      Neg_ras <- plot_ras[[Plot]]
      Neg_ras[which(values(Neg_ras) >= 0)] <- NA
      Pos_ras <- plot_ras[[Plot]]
      Pos_ras[which(values(Pos_ras) < 0)] <- NA
      if(Scaled == FALSE | Scaled == TRUE & plot == 1){
        mainTit <- Titles[Plot-1][[1]]
      }else{
        mainTit <- ""
      }
      plot(Neg_ras, col=col.signeg, colNA = "black", legend=FALSE, axes=FALSE, main = mainTit)
      if(Plot == 2 & min(values(plot_ras[[Plot]]), na.rm = TRUE) >= 0){
        plot(Pos_ras, col=col.sigpos, colNA = "black", legend=FALSE, axes=FALSE, add=TRUE)
      }else{
        plot(Pos_ras, col=col.sigpos, legend=FALSE, axes=FALSE, add=TRUE) 
      }
      min <- min(values(plot_ras[[Plot]]), na.rm = TRUE)
      max <- max(values(plot_ras[[Plot]]), na.rm = TRUE)
      range <- abs(min) + max
      if(length(Region) > 1 & Scaled == FALSE | Scaled == TRUE & plot == length(Region)){ # plot legend only for last raster stack if we scale
        minsegment <- .91 * abs(min)/range
        smallplotxpos <- c(.02+minsegment,.93,.08,.105) # where to put colour scales
        smallplotxneg <- c(0.02,.05+minsegment,.08,.105) # where to put colour scales
        plot(Neg_ras, legend.only=TRUE, col=col.signeg, colNA = "black", smallplot=smallplotxneg, horizontal = TRUE, axis.args=list(cex.axis=1.5))
        if(Plot == 2 & min(values(plot_ras[[Plot]]), na.rm = TRUE) >= 0){
          plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=c(.02, .93,.08,.105), horizontal = TRUE, axis.args=list(cex.axis=1.5))
        }else{
          plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=smallplotxpos, horizontal = TRUE, axis.args=list(cex.axis=1.5))
        }
      }else{
        minsegment <- .88 * abs(min)/range
        smallplotxpos <- c(.0+minsegment,.92,.15,.185) # where to put colour scales
        smallplotxneg <- c(0.05,.05+minsegment,.15,.185) # where to put colour scales
        if(Scaled == FALSE){
        plot(Neg_ras, legend.only=TRUE, col=col.signeg, colNA = "black", smallplot=smallplotxneg, horizontal = TRUE, axis.args=list(cex.axis=1.5))
        if(Plot == 2 & min(values(plot_ras[[Plot]]), na.rm = TRUE) >= 0){
          plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=1.5))
        }else{
          plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=smallplotxpos, horizontal = TRUE, axis.args=list(cex.axis=1.5))
        } 
        }
      }
    }
  }
  dev.off()
  
  ##------- VARIANCE PARTITIONING -------
  ## Plotting Setup
  # col.varpar1 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "wildfire")[30]
  # col.varpar2 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen")[50]
  # col.varpar3 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "jon_snow")[1]
  # col.list <- list(col.varpar1, col.varpar2, col.varpar3) 
  col.list <- as.list(viridis(1000)[c(1, 500, 1000)])
  
  ## Data
  MaxVar <- rep(NA, length(Region))
  for(plot in 1:length(Region)){
    Alter_ras <- Rasters[[plot]][[7:10]]
    Alter_ras[2] <- 0
    values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
    cells <- order(values(Alter_ras[[1]]))
    `%nin%` = Negate(`%in%`) # create a 'not in' statement
    cells <- cells[which(cells %nin% which(values(Alter_ras[[1]])>quantile(values(Alter_ras[[1]]), .95, na.rm = TRUE)))]
    MaxVar[plot] <- max(Alter_ras[cells], na.rm = TRUE)
  }
  VarPars <- list()
  for(plot in 1:length(Region)){
    Alter_ras <- Rasters[[plot]][[7:10]]
    Alter_ras[2] <- 0
    values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
    cells <- order(values(Alter_ras[[1]]))
    `%nin%` = Negate(`%in%`) # create a 'not in' statement
    cells <- cells[which(cells %nin% which(values(Alter_ras[[1]])>quantile(values(Alter_ras[[1]]), .95, na.rm = TRUE)))]
    plot_df <- data.frame(Data =NA, Cell =NA, Variance =NA)
    Idents <- c("Total", "t-1", "Shared", "Qsoil")
    for(i in 1:4){
      if(i > 1){
        plot_df1 <- data.frame(Data = values(Alter_ras[[i]])[cells],
                               Cell = 1:length(cells), 
                               Variance = rep(Idents[i], length(cells)))
        plot_df <- rbind(plot_df, plot_df1)
      }
    }
    Lims <- c(0, max(plot_df$Data, na.rm = TRUE))
    plot_df <- na.omit(plot_df)
    if(length(Region) == 1){
      p <- ggplot(data = plot_df, aes(y = Data, x = Cell, fill = Variance)) + geom_bar(stat = "identity") + theme_bw(base_size= 25) + xlab("Raster Cells") + ylab("Variance") + scale_fill_manual(values=c(col.list[[3]], col.list[[2]], col.list[[1]])) + ylim(Lims)
    }else{
      if(Scaled == TRUE){
        p <- ggplot(data = plot_df, aes(y = Data, x = Cell, fill = Variance)) + geom_bar(stat = "identity") + theme_bw(base_size= 25) + scale_fill_manual(values=c(col.list[[3]], col.list[[2]], col.list[[1]])) + ylim(Lims) + xlab("Raster Cells") + ylab("Variance") + ggtitle(Region[plot]) + scale_y_continuous(limits = c(0, max(MaxVar)))
      }else{
        p <- ggplot(data = plot_df, aes(y = Data, x = Cell, fill = Variance)) + geom_bar(stat = "identity") + theme_bw(base_size= 25) + scale_fill_manual(values=c(col.list[[3]], col.list[[2]], col.list[[1]])) + ylim(Lims) + xlab("Raster Cells") + ylab("Variance") + ggtitle(Region[plot])
      }
      
    }
    VarPars[[plot]] <- p
  }
  ## Saving Files
  if(length(Region) ==1){
    jpeg(file=paste(Dir.Plots, "/", paste(Region,collapse=""), "_VarPar", SoilLayer, "_", Scaled, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 100)
    print(VarPars[[1]])
    dev.off()
  }else{
    height <- 22*round(length(Region)/2)
    width <- 42
    jpeg(file=paste(Dir.Plots, "/", paste(Region,collapse=""), "_VarPar", SoilLayer, "_", Scaled, ".jpeg", sep = ""), width = width, height = height, units = "cm", quality = 100, res = 100)
    
    if(length(Region)==2){
      p <- ggarrange(VarPars[[1]], VarPars[[2]], ncol=2, 
                common.legend = TRUE, legend="bottom",
                labels = "AUTO")
      print(p)
      
    }
    if(length(Region)==3){
      leg <- VarPars[[1]] + theme_bw(base_size = 15) + guides(fill = guide_legend(override.aes = list(size = 15)))
      legend <- cowplot::get_legend(leg)
      grid.arrange(VarPars[[1]]+ theme(legend.position = "none"),
                   VarPars[[2]]+ theme(legend.position = "none"), 
                   VarPars[[3]]+ theme(legend.position = "none"),
                   legend, 
                   ncol=2,
                   left = textGrob("Explained Variance", rot = 90, vjust = 1),
                   bottom = textGrob("Raster Cells", rot = 0, vjust = 1))
    }
    if(length(Region)==4){
      p <- ggarrange(VarPars[[1]], VarPars[[2]], VarPars[[3]], VarPars[[4]], ncol=2, nrow = 2, 
                common.legend = TRUE, legend="bottom",
                labels = "AUTO")
      print(p)
    }
    if(length(Region)==5){
      leg <- VarPars[[1]] + theme_bw(base_size = 15) + guides(fill = guide_legend(override.aes = list(size = 15)))
      legend <- cowplot::get_legend(leg)
      grid.arrange(VarPars[[1]]+ theme(legend.position = "none"),
                   VarPars[[2]]+ theme(legend.position = "none"), 
                   VarPars[[3]]+ theme(legend.position = "none"),
                   VarPars[[4]]+ theme(legend.position = "none"),
                   VarPars[[5]]+ theme(legend.position = "none"),
                   legend, 
                   ncol=2,
                   left = textGrob("Explained Variance", rot = 90, vjust = 1),
                   bottom = textGrob("Raster Cells", rot = 0, vjust = 1))
    }
    dev.off()
  }
  
  # Maps
  col.varpar1 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "wildfire")
  col.varpar2 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen")
  col.varpar3 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "jon_snow")
  col.list <- list(col.varpar1, col.varpar2, col.varpar3)
  
  for(plot in 1:length(Region)){
  jpeg(file=paste(Dir.Plots, "/", Region[plot], "_VarParMap", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(2,2))
  Alter_ras <- Rasters[[plot]][[7:10]]
  Alter_ras[2] <- 0
  values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
  for(i in 2:4){
    col.varpar <- col.list[[i-1]]
    plot(Alter_ras[[i]], col=col.varpar, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[i-1])
    plot(Alter_ras[[i]], legend.only=TRUE, col=col.varpar, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  }
  dev.off()
  }
} # end of Fun_Plot