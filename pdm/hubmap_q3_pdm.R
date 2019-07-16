library(imager)
library(colorspace)
library(ggplot2)
library(Seurat)
library(ggforce)

# load the lung image
lung <- load.image("~/xfer/lung.png")
lung <- R(lung)

lung.gg <- as.data.frame(lung)

# binarize to 0-1
lung.gg[, "value"][lung.gg[, "value"] > 0.5] <- 1
lung.gg[, "value"][lung.gg[, "value"] < 0.5] <- 0
lung.gg[, "value"] <- (lung.gg[, "value"] - 1) * -1

lung.binary <- as.cimg(lung.gg)
# compute distance transformation
lung.dist <- as.data.frame(distance_transform(lung.binary, 1))

# make background grey
lung.dist[lung.dist$value == 0, "value"] <- NA
lung.dist <- subset(lung.dist, x < 750 & x > 100)
lung.dist$v2 <- lung.dist$value / max(lung.dist$value, na.rm = T)

p1 <- ggplot(lung.dist, aes(x,y)) + geom_raster(aes(fill = value)) +
  ylim(c(min(lung.dist$y), max(lung.dist$y))) + xlim(c(min(lung.dist$x), max(lung.dist$x))) + 
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())+ 
  scale_fill_continuous_sequential(palette = "Purples 3", end = 1, begin = 0.4)

origin_1 <- c(390, 607)
p1 <- p1 + geom_point(data = data.frame(x = origin_1[1], y = origin_1[2]), aes(x = x, y = y), color = "firebrick", size = 3)

# landmarkx_1 <- c(575, 636)
landmarkx_1 <- c(716, 652)
landmarky_1 <- c(485, 125)

vx_1 <- landmarkx_1 - origin_1
vy_1 <- landmarky_1 - origin_1

point_distance <- function(par, point, origin, vx, vy) {
  p1 <- origin+vx*(par[1])+vy*(par[2])
  return(sqrt(sum((p1-point)^2)))
}

lung2 <- load.image("~/xfer/lung3.png")
lung2 <- R(lung2)
lung2.gg <- as.data.frame(lung2)

# binarize to 0-1
lung2.gg[, "value"][lung2.gg[, "value"] > 0.5] <- 1
lung2.gg[, "value"][lung2.gg[, "value"] < 0.5] <- 0
lung2.gg[, "value"] <- (lung2.gg[, "value"] - 1) * -1

lung2.binary <- as.cimg(lung2.gg)
# compute distance transformation
lung2.dist <- as.data.frame(distance_transform(lung2.binary, 1))

# make background grey
lung2.dist[lung2.dist$value == 0, "value"] <- NA
lung2.dist <- subset(lung2.dist,x<500&x>10)
lung2.dist$v2 <- lung2.dist$value / max(lung2.dist$value,na.rm = T)

# plot 
p2 <- ggplot(lung2.dist, aes(x,y)) + geom_raster(aes(fill = value)) +
  ylim(c(min(lung2.dist$y), max(lung2.dist$y))) + xlim(c(min(lung2.dist$x), max(lung2.dist$x))) + 
  scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),trans=scales::reverse_trans())+ 
  scale_fill_continuous_sequential(palette = "Purples 3", end = 1, begin = 0.4)

origin_2 <- c(237,553)
p2 <- p2 + geom_point(data = data.frame(x = origin_2[1], y = origin_2[2]), aes(x = x, y = y), color = "firebrick", size = 3)

landmarkx_2 <- c(476,621)
landmarky_2 <- c(355,165)

vx_2 <- landmarkx_2 - origin_2
vy_2 <- landmarky_2 - origin_2


pdm_distance <- function(par,v,origin,distmap,targetpdm) {
  newpoint <- round(origin_2 + v*par[1])
  newpdm <- subset(distmap, x==newpoint[1] & y==newpoint[[2]])$v2
  return(abs(targetpdm-newpdm))
}

steps <- 100
circle.df <- data.frame(
  x0 = rep(0, steps),
  y0 = rep(0, steps),
  r = (steps:1)/steps,
  fill = (steps:1)/steps
)

### randomly generate points to transfer
set.seed(42)
site1 <- sample(x = rownames(na.omit(lung.dist)), size = 1)
plot.list <- list()
i = 1

for(s in site1) {
  print(i)
  s1 <- as.numeric(lung.dist[s, 1:2])
  pdm1 <- subset(lung.dist,x==s1[1] & y==s1[[2]])$v2
  point1 <- s1 - origin_1
  proj_1 <- as.numeric((point1%*%vx_1)/(vx_1%*%vx_1))*vx_1 + origin_1
  result <- optim(par=c(-1,1), fn = point_distance, point = s1, origin=origin_1,vx= vx_1, vy = vy_1)
  site1projection <- origin_2+vx_2*(result$par[1])+vy_2*(result$par[2])
  sitevec_2 <- site1projection - origin_2
  result2 <- try(expr = optim(par=c(0.8), fn = pdm_distance, v = sitevec_2, origin = origin_2, distmap = lung2.dist,targetpdm = pdm1), silent = T)
  if(class(x = result2) == "try-error") next
  site2 <- round(origin_2 + sitevec_2*result2$par[1])
  p3 <- p1 + annotate("point", x = s1[1], y = s1[2], colour = "darkorange1", size = 3)
  p4 <- p2 + annotate("point", x = site2[1], y = site2[2], colour = "darkorange1", size = 3)
  
  xpos <- (s1 - origin_1)[1]
  ypos <- -(s1 - origin_1)[2]
  
  pos_vector <- matrix(data = c(xpos, ypos), ncol = 1)
  vec_len <- sqrt(crossprod(pos_vector))
  
  norm_vec <- (pos_vector /  as.numeric(vec_len)) * (1 - pdm1)
  
  df <- as.data.frame(t(norm_vec))
  
  plot_pdm <- ggplot(df, aes(x = 0, y = 0)) +
    geom_circle(data = circle.df, aes(x0 = x0, y0 = y0, r = r, fill = fill, color = fill), inherit.aes = FALSE) +
    #geom_segment(aes(xend = V1, yend = V2), arrow = arrow(length = unit(0.3, "cm")))+
    xlim(c(-1, 1)) +
    ylim(c(-1, 1)) +
    geom_point(data = data.frame(x = 0, y = 0), aes(x = x, y = y), color = "firebrick", size = 3) + 
    coord_fixed() + 
    scale_fill_continuous_sequential(palette = "Purples 3", end = 0.4, begin = 1) +
    scale_color_continuous_sequential(palette = "Purples 3", end = 0.4, begin = 1) +
    NoLegend() + NoAxes()
  
  c1 <- plot_grid(plotlist = list(p3 + NoLegend(), plot_pdm, p4 + NoLegend()), ncol = 3)
  plot.list[[i]] <- c1
  i <- i + 1
}
