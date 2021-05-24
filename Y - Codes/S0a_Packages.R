install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c("automap", "doParallel", "foreach", "gameofthrones", "ggplot2", "gimms", "mapview", 
                 "gridExtra", "ncdf4", "pracma", "raster", "rgbif", "rgdal", "sp", "vegan", "RStoolbox", "forcats", "grid", "rasterVis", "cowplot", "BIEN", "ggpubr", "viridis", "modEvA", "nlme", "car", "jcolors", "KrigR", "forecast")

`%nin%` = Negate(`%in%`)
if("KrigR" %nin% installed.packages()){
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  devtools::install_github("https://github.com/ErikKusch/KrigR") 
}

sapply(package_vec, install.load.package)
