################################################################################
####--------------- Fun_Plot [Region, Scaled] # all things vegetation memory (maps, varpar, system dynamics)
Fun_Plot <- function(Region, Scaled = FALSE){
  ####--------------- FILE SELECTION ----------------
  Raster <- brick(paste(Dir.Memory, "/", Region, ".nc", sep=""))
  plot(Raster[[15]])
  Back <- raster(paste(Dir.KrigCov, "Co-variates_NativeResolution.nc", sep="/"),
                 varname = "Elevation")
  Countries <- readOGR(Dir.Mask, "ne_50m_admin_0_countries", verbose = FALSE)
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
      values(Raster[[rasterscale]])[which(!is.na(values(Raster[[rasterscale]])))[101]] <- max(maxs)
      values(Raster[[rasterscale]])[which(!is.na(values(Raster[[rasterscale]])))[100]] <- min(mins)
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
    par(mfrow = c(3,2), mai = c(.9, 0.2, 0.2, 0.2))
    smaplot <- c(.05, .95, .27, .3)
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
          smallplotxpos <- c(.05+minsegment,.95,.27,.3) # where to put colour scales
          smallplotxneg <- c(0.05,.05+minsegment,.27,.3) # where to put colour scales
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
  P_Box <- ggplot(data = plot_df, aes(y = Data, x = Variance, fill = Variance)) + geom_boxplot() + theme_bw(base_size= 25) + xlab("Raster Cells") + ylab("Variance") + 
    scale_fill_manual(values=c(col.list[[1]], col.list[[2]], col.list[[3]], col.list[[4]],
                               col.list[[5]], col.list[[6]], col.list[[7]])) + 
    theme(axis.text.x = element_text(size=7, angle=-30))
  P_VarsA <- P_Diag +
    annotation_custom(
      grob = ggplotGrob(P_Box + theme(legend.position = "none", 
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank())),
      xmin = 0,
      xmax = ceiling(quantile(plot_df$Cell, .7)),
      ymin = .15,
      ymax = .5
    )
  P_VarsB <- P_Box +
    annotation_custom(
      grob = ggplotGrob(P_Diag + theme(legend.position = "none", 
                                       axis.title.x = element_blank(),
                                       axis.title.y = element_blank())),
      xmin = 1.1,
      xmax = 6.9,
      ymin = .1,
      ymax = .25
    )
  ## Saving Plots
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParDiag.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 100)
  P_Diag
  dev.off()
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParBox.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 100)
  P_Box
  dev.off()
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParA.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 100)
  P_VarsA
  dev.off()
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParB.jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 100)
  P_VarsB
  dev.off()
  
  # Maps
  col.varpar1 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "wildfire")
  col.varpar2 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen")
  col.varpar3 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "jon_snow")
  col.varpar4 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "margaery")
  col.list <- list(col.varpar1, col.varpar1, col.varpar2, col.varpar3,
                   col.varpar4, col.varpar4, col.varpar4, col.varpar4)
  Titles <- c("Total Explained Variance", "NDVI[t-1]", "Qsoil1", "Tair", "NDVI[t-1] + Qsoil1",
              "NDVI[t-1] + Tair", "Qsoil1 + Tair", "Shared by all")
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParMap.jpeg", sep = ""), width = 22, height = 32, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(4,2))
  Alter_ras <- Raster[[7:14]]
  Alter_ras[2] <- 0
  values(Alter_ras)[which(values(Alter_ras) < 0)] <- 0
  for(i in 1:8){
    col.varpar <- col.list[[i]]
    plot(Alter_ras[[i]], col=col.varpar, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[i])
    plot(Alter_ras[[i]], legend.only=TRUE, col=col.varpar, smallplot=c(.05, .95, .27, .3), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  }
  dev.off()
} # end of Fun_Plot
