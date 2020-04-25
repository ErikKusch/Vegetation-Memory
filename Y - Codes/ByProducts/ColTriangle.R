# Coordinates of the triangle
tri <- rbind(sin(0:2*2/3*pi), cos(0:2*2/3*pi))

# Function for calculating the color of a set of points `pt`
# in relation to the triangle
tricol <- function(pt, sharpness=2){
  require(splancs)
  RGB <- sapply(1:3, function(i){
    a <- sweep(pt, 2, tri[,i])
    b <- apply(tri[,-i], 1, mean) - tri[,i]
    sharpness*((a %*% b) / sum(b^2))-sharpness+1
  })
  RGB[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
  do.call(rgb, unname(as.data.frame(pmin(pmax(RGB, 0), 1))))
}

# Plot
res <- 1000                         # Resolution
xi <- seq(-1, 1, length=res)        # Axis points
yi <- seq(-.8, 1.2, length=res)
x <- xi[1] + cumsum(diff(xi))       # Midpoints between axis points
y <- yi[1] + cumsum(diff(yi))
xy <- matrix(1:(length(x)*length(y)), length(x))
image(xi, yi, xy, col=tricol(as.matrix(expand.grid(x,y))), useRaster=TRUE)
lines(tri[1,c(1:3,1)], tri[2,c(1:3,1)], type="l")
