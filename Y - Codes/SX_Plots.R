################################################################################
####--------------- Fun_Plot [Region, Scaled] # all things vegetation memory (maps, varpar, system dynamics)
Fun_Plot <- function(Region, Parameter = NULL){
  ####--------------- FILE SELECTION ----------------
  Data <- read.csv(file.path(Dir.Memory, paste0(Region, ".csv")))
  
  ####--------------- PLOTTING ----------------
  if(Parameter == "Memory"){
    Data2 <- gather(Data[,c(4:8, 19)], Parameter, Value, NDVI:Lag.Tair, factor_key = TRUE)
    ggplot(data = Data2, aes(x = Parameter, y = Value)) + 
      geom_boxplot() + 
      facet_wrap(~Detrending, ncol = 3) + 
      theme_bw() + labs(x = "Memory Component", y = "")
    ggsave(filename = file.path(Dir.Plots, paste0("Detrending_", Parameter, ".png")), width = 16, height = 9)
  }
  if(Parameter == "Variances"){
    Data2 <- gather(Data[,c(9:16, 19)], Parameter, Value, V_NDVI:V_Shared, factor_key = TRUE)
    ggplot(data = Data2, aes(x = Parameter, y = Value)) + 
      geom_boxplot() + 
      facet_wrap(~Detrending, ncol = 3) + 
      theme_bw() + labs(x = "Memory Component", y = "Variance Explained")
    ggsave(filename = file.path(Dir.Plots, paste0("Detrending_", Parameter, ".png")), width = 16, height = 9)
  }
  if(Parameter != "Memory" & Parameter != "Variances"){
    ggplot(data = Data, aes(x = Detrending, y = Data[,Parameter])) + 
      geom_boxplot() + 
      theme_bw() + labs(y = Parameter)
    ggsave(filename = file.path(Dir.Plots, paste0("Detrending_", Parameter, ".png")), width = 16, height = 9)
  }
  
} # end of Fun_Plot