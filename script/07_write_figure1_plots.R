# Write maps used in the illustration of Fig 01.
#...............................................
#...............................................
library("raster")
library("sf")
library("smoothr")
library("dismo")
library("ggplot2")

#...............................................
#...............................................
# define the extent 
#e <- extent(-86, -85.3, 12.8, 13.42)
e <- extent(-85.95, -85.3, 12.8, 13.41)

#...............................................
#...............................................
# this file is an illustration for the idea of centralized approach,
# recommendations based on administrative regions
l <- getData("GADM", country = "NIC", level=1)

# crop file for the extent
l <- crop(l, e)

l <- st_as_sf(l)

l$NAME_1 <- c("Var A", "Var A", "Var B")

colors <- c("#253dac","#70e072")

p <- ggplot(l) +
  geom_sf(aes(fill = as.factor(NAME_1)),
          lwd = 0) + 
  scale_fill_manual(values = colors,
                    labels = c("Var A", "Var B")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 15, colour="black"))

p

ggsave("output/fig1/centralized.svg", 
       plot = p,
       width = 11, 
       height = 11,
       units = "cm")

#...............................................
#...............................................
# this file is an illustration for the idea of decentralized approach
# recommendations based on environmental conditions
d <- raster("data/fig1.tif")
plot(d)
plot(e, add = TRUE)
d <- crop(d, e)

plot(d)


# replace some values
extf <- extent(c(-85.95,-85.65,12.8,13.14))
extf2 <- extent(c(-85.8,-85.42,13.25,13.41))
extf3 <- extent(c(-85.5,-85.40,13.13,13.25))
extf4 <- extent(c(-85.7,-85.41,12.8,12.9))
extf5 <- extent(c(-85.55,-85.38,13,13.12))
extf6 <- extent(c(-85.8,-85.68,12.8,12.88))
extf7 <- extent(c(-85.5,-85.43,13.22,13.3))

plot(d)
#plot(extf7, add = TRUE)
plot(extf6, add = TRUE)
plot(extf5, add = TRUE, col = "red")
plot(extf, add = TRUE, col = "grey")
plot(extf2, add = TRUE, col = "black")
plot(extf3, add = TRUE, col = "pink")
plot(extf4, add = TRUE, col = "green")

# Change values
d[extf][d[extf] == 3] <- 4
d[extf][d[extf] == 1] <- 4

d[extf2][d[extf2] == 2] <- 1

d[extf3][d[extf3] == 1] <- 2
d[extf3][d[extf3] == 3] <- 2

d[extf4][d[extf4] == 1] <- 2
d[extf4][d[extf4] == 3] <- 2

d[extf5][d[extf5] == 3] <- 2

d[extf6][d[extf6] == 2] <- 4


d4 <- as.data.frame(d, xy = TRUE)
d4 <- d4[d4$fig1==4,]
d4 <- d4[d4$y > 12.85,]
d4 <- d4[d4$x > -85.84,]
d4 <- dismo::convHull(d4)
d4 <- d4@polygons
r <- raster(d)
d4 <- rasterize(d4, r)

plot(d)
plot(d4, add = TRUE)


d[d4][d[d4] !=10] <- 4

plot(d)

s <- rasterToPolygons(d, dissolve = TRUE)
s <- st_as_sf(s)
names(s)[1] <- "value"
s <- smooth(s, method = "chaikin", refinements = 1)
s <- smooth(s, method = "ksmooth", smoothness = 1)
plot(s)

colors <- c("#253dac",
            "#70e072",
            "#f1c73e",
            "#2a4d0c")
#colors <- sample(colors, 4)

p2 <- 
  ggplot(s) + 
  geom_sf(aes(fill = factor(value, levels = c("1","4","3","2"))),
          lwd = 0) +
  scale_fill_manual(values= colors,
                    name = NULL,
                    labels = c("Var A", "Var B", "Var C", "Var D")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 15, colour="black"))

p2

ggsave("output/fig1/decentralized.svg", 
       plot = p2,
       width = 11,
       height = 11, 
       units = "cm")

# remove the layer file
# this can be downloaded again from GADM
file.remove("gadm36_NIC_1_sp.rds")
