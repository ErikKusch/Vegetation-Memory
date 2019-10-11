################################################################################
rm(list = ls()) # clearing environment
####--------------- PACKAGES ----
source("Y - Codes/S0a_Packages.R") # loading packages
####--------------- DIRECTORIES ----
source("Y - Codes/S0b_Directories.R") # setting directories
####--------------- FUNCTIONS ----
source("Y - Codes/S0c_Functions.R") # Loading miscellaneous functions
####--------------- VARIABLE VECTORS ----
ClimVars = list("Qsoil1_mean", "Qsoil2_mean", "Qsoil3_mean", "Qsoil4_mean")
ClimVars2 = list("Tair_mean", "Tair_mean", "Tair_mean", "Tair_mean")

################################################################################
####--------------- PLOTTING FUNCTION ----
Fun_Plot <- function(Region, SoilLayer){
  ####--------------- FILE SELECTION ----------------
  Files <- list.files(Dir.Memory)[grep(pattern = Region, list.files(Dir.Memory))]
  File <- Files[grep(pattern = paste("Qsoil", SoilLayer, sep=""), Files)]
  Memory <- brick(paste(Dir.Memory,File, sep="/"))
  ####--------------- MISC ----------------
  SR_Titles <- list("NDVI[t-1]", "Air Temperature", "Soil Moisture (0-7cm)", "Soil Moisture (7-28cm)", "Soil Moisture (28-100cm)", "Soil Moisture (100-255cm)")
  Titles <- c(paste(SR_Titles[1], "(Intrinsic Memory)"), paste(SR_Titles[2+SoilLayer], "(Inverse Resistance)"), paste(SR_Titles[2], "(Inverse Resistance)"))
  TitlesShort <- c(SR_Titles[1][[1]], SR_Titles[2+SoilLayer][[1]], SR_Titles[2][[1]])
  ####--------------- PLOTTING ----------------
  ##------- MEMORY COMPONENTS -------
  ## Plotting Setup
  col.signeg <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "targaryen2")
  col.sigpos <- got(n = 100, alpha = 1, begin = 0, end = 1, direction = -1, option = "tyrell")
  col.nonsig <- colorRampPalette(c("grey"))(1)
  col.lags <- got(n = 12, alpha = 1, begin = 0, end = 1, direction = 1, option = "daenerys")
  ## Data
  plot_ras <- Memory
  ## Saving Files
  # TriColour
  Triplot <- abs(plot_ras[[4:6]])
  TriMap <- ggRGB(img = Triplot, r = 3, g = 1, b = 2, scale = max(maxValue(Triplot)), stretch = 'none') + ggtitle(label = "Relative Importance of Vegetation Responses") + ylab("Latitude [°]") + xlab("Longitude [°]") + theme_bw(base_size = 15)
  # + scale_colour_manual(name = '', values =c('black'='black','red'='red'), labels = c('c2','c1'))
  ggframe <- data.frame(Values = as.vector(Triplot),
             Idents = rep(TitlesShort, each = dim(values(Triplot))[1]),
             Order = rep(1:3, each = dim(values(Triplot))[1]))
  ggframeSmall <- na.omit(ggframe)
  TriBars <- ggplotGrob(ggplot(ggframeSmall, aes(x = fct_reorder(Idents,Order), y = Values, fill = Idents)) + geom_boxplot() + theme_minimal(base_size = 15) + guides(fill=guide_legend(title="")) + theme(legend.position = c(.6, .9)) + ylab("Coefficients") + 
    theme(panel.background = element_rect(fill = "lightgrey")) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    )
  GGCombined <- TriMap +
    annotation_custom(
      grob = TriBars,
      xmin = -0.5, # these coordinates need altering for different regions
      xmax = 5,
      ymin = 34.6,
      ymax = 38.75
    )
  ggsave(file=paste(Dir.Plots, "/", Region, "_RelImportance", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
  # Two-by-Two
  jpeg(file=paste(Dir.Plots, "/", Region, "_Model", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(2,2))
  plot(plot_ras[[1]], col = col.lags, colNA = "black", legend=FALSE, axes=FALSE, main = paste(SR_Titles[2+SoilLayer][[1]], "(Memory Length)"))
  plot(plot_ras[[1]], legend.only=TRUE, col=col.lags, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  for(Plot in 4:6){
    Neg_ras <- plot_ras[[Plot]]
    Neg_ras[which(values(Neg_ras) >= 0)] <- NA
    Pos_ras <- plot_ras[[Plot]]
    Pos_ras[which(values(Pos_ras) < 0)] <- NA
    plot(Neg_ras, col=col.signeg, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[Plot-3][[1]])
    if(Plot == 4 & min(values(plot_ras[[Plot]]), na.rm = TRUE) >= 0){
      plot(Pos_ras, col=col.sigpos, colNA = "black", legend=FALSE, axes=FALSE, add=TRUE)
    }else{
      plot(Pos_ras, col=col.sigpos, legend=FALSE, axes=FALSE, add=TRUE) 
    }
    min <- min(values(plot_ras[[Plot]]), na.rm = TRUE)
    max <- max(values(plot_ras[[Plot]]), na.rm = TRUE)
    range <- abs(min) + max
    minsegment <- .88 * abs(min)/range
    smallplotxpos <- c(.05+minsegment,.93,.15,.185) # where to put colour scales
    smallplotxneg <- c(0.05,.05+minsegment,.15,.185) # where to put colour scales
    plot(Neg_ras, legend.only=TRUE, col=col.signeg, colNA = "black", smallplot=smallplotxneg, horizontal = TRUE, axis.args=list(cex.axis=2.5))
    if(Plot == 4 & min(values(plot_ras[[Plot]]), na.rm = TRUE) >= 0){
      plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=2.5))
    }else{
      plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=smallplotxpos, horizontal = TRUE, axis.args=list(cex.axis=2.5))
    }
  }
  dev.off()
  
  ##------- VARIANCE PARTITIONING -------
  ## Plotting Setup
  col.varpar1 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "wildfire")
  col.varpar2 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "targaryen")
  col.varpar3 <- got(n = 100, alpha = 1, begin = 0.2, end = 1, direction = -1, option = "jon_snow")
  col.list <- list(col.varpar1, col.varpar2, col.varpar3) 
  ## Data
  Alter_ras <- Memory[[7:10]]
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
  ## Saving Files
  # Stacked Bars
  p <- ggplot(data = plot_df, aes(y = Data, x = Cell, fill = Variance)) + geom_bar(stat = "identity") + theme_bw(base_size= 45) + xlab("Raster Cells") + ylab("Variance") + scale_fill_manual(values=c(col.list[[3]][1], col.list[[2]][50], col.list[[1]][30])) + ylim(Lims)
  ggsave(file=paste(Dir.Plots, "/", Region, "_VarPar", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
  # Maps
  jpeg(file=paste(Dir.Plots, "/", Region, "_VarParMap", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  par(mfrow=c(2,2))
  for(i in 2:4){
    col.varpar <- col.list[[i-1]]
    plot(Alter_ras[[i]], col=col.varpar, colNA = "black", legend=FALSE, axes=FALSE, main = Titles[i-1])
    plot(Alter_ras[[i]], legend.only=TRUE, col=col.varpar, smallplot=c(.05, .93, .15, .185), horizontal = TRUE, axis.args=list(cex.axis=2.5))
  } 
  dev.off()
  
  ##------- SYSTEM DYNAMICS -------
  ## Preparation
  # assuming a negative, linear relationship between recovery rates and intrinsic vegetation memory, the higher this proxy, the faster the recovery at constant intrinsic vegetation memory
  InverseProx <- 10
  Equilibrium <- 0
  plot_ras <- Memory[[c(1,4:5)]]
  # Plot Raster for Selection
  Neg_ras <- plot_ras[[3]]
  Neg_ras[which(values(Neg_ras) >= 0)] <- NA
  Pos_ras <- plot_ras[[3]]
  Pos_ras[which(values(Pos_ras) < 0)] <- NA
  plot(Neg_ras, col=col.signeg, colNA = "black", legend=FALSE, axes=FALSE, main = "Select Cells By Clicking and Finish with ESCAPE. Select entire region without clicking on raster and hitting ESCAPE.")
  plot(Pos_ras, col=col.sigpos, legend=FALSE, axes=FALSE, add=TRUE) 
  min <- min(values(plot_ras[[3]]), na.rm = TRUE)
  max <- max(values(plot_ras[[3]]), na.rm = TRUE)
  range <- abs(min) + max
  minsegment <- .88 * abs(min)/range
  smallplotxpos <- c(.05+minsegment,.93,.15,.185) # where to put colour scales
  smallplotxneg <- c(0.05,.05+minsegment,.15,.185) # where to put colour scales
  plot(Neg_ras, legend.only=TRUE, col=col.signeg, colNA = "black", smallplot=smallplotxneg, horizontal = TRUE, axis.args=list(cex.axis=2.5))
  plot(Pos_ras, legend.only=TRUE, col=col.sigpos, smallplot=smallplotxpos, horizontal = TRUE, axis.args=list(cex.axis=2.5))
  ## Data
  # Prompt Cell Selection
  Click_Pos <- click(plot_ras[[3]], cell = TRUE)$cell
  if(is.null(Click_Pos)){
    Data_Pos <- which(!is.na(values(plot_ras[[1]])))
  }else{
    Data_Pos <- adjacent(plot_ras[[3]], cells=Click_Pos, directions=16, pairs=FALSE, include = TRUE)  
  }
  # Prepare Empty Data Frame and Progress Bar
  plot_df <- data.frame(State = NA, Time = NA, Response = NA)
  counter <- 0
  print("Extracting Data")
  pb <- txtProgressBar(min = 0, max = length(Data_Pos), style = 3)
  for(i in Data_Pos){
    if(is.na(plot_ras[[3]][i])){#coastline adjacent cells
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
      next()
      }
    if(plot_ras[[3]][i] < 0){
      Identifier <- "Attenuation"
    }else{
      Identifier <- "Resonance"
    }
    
    States <- c(Equilibrium, plot_ras[[3]][i], 
                Equilibrium)
    Times <- c(2,2+plot_ras[[1]][i], 
               2+plot_ras[[1]][i]+abs(plot_ras[[3]][i]*plot_ras[[2]][i]/InverseProx))
    Response <- rep(Identifier, 3)
    
    plot_df <- rbind(plot_df, 
                     data.frame(State = States, Time = Times, Response = Response))
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
  plot_df <- na.omit(plot_df)
  plot_df <- rbind(plot_df, data.frame(State=rep(Equilibrium, 2),
                              Time =c(0,2),
                              Response = rep("Equilibrium",2)) )
  ## Plotting, these one needs adding of vertical lines and text!!!
  # dynamics
  Lines_gg <- ggplot(plot_df, aes(x = Time, y = State, col = Response)) + 
    geom_line(data = plot_df[which(plot_df$Response == "Equilibrium"),]) + 
    stat_smooth(level = 0.95, fullrange = FALSE, data = plot_df[which(plot_df$Response != "Equilibrium"),]) + 
    theme_bw(base_size = 15) + scale_color_manual(values=c("red", "black", "forestgreen")) + 
    geom_vline(aes(xintercept = 2), col = "blue") + 
    geom_text(aes(x=1.8, label="Positive Soil Moisture Anomaly", y=1.5), colour="blue", angle=90, text=element_text(size=11))
# insert map
  Mini_ras <- plot_ras[[3]]
  Mini_ras[Data_Pos[which(plot_ras[[3]][Data_Pos] < 0)]] <- 0
  Mini_ras[Data_Pos[which(plot_ras[[3]][Data_Pos] >= 0)]] <- 10000
  Mini_ras[which(values(Mini_ras) != 10000 & values(Mini_ras) != 0)] <- -10000
  Miniras_gg <- gplot(Mini_ras) + geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = 'grey', mid = "red", high = 'forestgreen', na.value = NA) +
    coord_equal() + theme_void(base_size = 15) + theme(legend.position = "none")
  # combining
  Combined_gg <- ggdraw() + draw_plot(Lines_gg) + draw_plot(Miniras_gg, 0.53, 0.65, 0.35, 0.3)
  ggsave(file=paste(Dir.Plots, "/", Region, "_Dynamics", SoilLayer, ".jpeg", sep = ""), width = 32, height = 22, units = "cm", quality = 100)
} # end of Fun_Plot

################################################################################
Fun_Plot(Region = "SWEurope", SoilLayer = 1)
Fun_Plot(Region = "SWEurope", SoilLayer = 2)