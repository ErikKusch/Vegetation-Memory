install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c("automap", "doParallel", "foreach", "gameofthrones", "ggplot2", "gimms", 
                 "gridExtra", "ncdf4", "pracma", "raster", "rgbif", "rgdal", "sp", "vegan", "xlsx", "RStoolbox", "forcats", "grid", "rasterVis", "cowplot", "BIEN", "ggpubr")
sapply(package_vec, install.load.package)
