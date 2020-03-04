install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c("automap", "doParallel", "foreach", "gameofthrones", "ggplot2", "gimms", "mapview", 
                 "gridExtra", "ncdf4", "pracma", "raster", "rgbif", "rgdal", "sp", "vegan", "RStoolbox", "forcats", "grid", "rasterVis", "cowplot", "BIEN", "ggpubr", "viridis", "modEvA")
sapply(package_vec, install.load.package)
