library(raster)
list.files()
setwd(mainDir)
D1 <- raster(list.files()[grep(pattern = "D01", list.files())])
D2 <- raster(list.files()[grep(pattern = "D2", list.files())])
D3 <- raster(list.files()[grep(pattern = "D3", list.files())])
D4 <- raster(list.files()[grep(pattern = "D4", list.files())])
D5 <- raster(list.files()[grep(pattern = "D5", list.files())])
D6 <- raster(list.files()[grep(pattern = "D6", list.files())])
D7 <- raster(list.files()[grep(pattern = "D7", list.files())])
D8 <- raster(list.files()[grep(pattern = "D8", list.files())])
D9 <- raster(list.files()[grep(pattern = "D9", list.files())])
D10 <- raster(list.files()[grep(pattern = "D10", list.files())])
D11 <- raster(list.files()[grep(pattern = "D11", list.files())])
D12 <- raster(list.files()[grep(pattern = "D12", list.files())])
D13 <- raster(list.files()[grep(pattern = "D13", list.files())])

Ds <- merge(D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13)
setwd(Dir.Gimms.Monthly)
Back <- brick(list.files(Dir.Gimms.Monthly)[1])[[6]]
setwd(mainDir)
BackC <- crop(Back, extent(Ds))

plot(BackC, col = "grey", legend = FALSE)
plot(Ds, add = TRUE)


