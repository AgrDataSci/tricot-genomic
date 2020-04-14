# Write maps used in the illustration of Fig 01.
#...............................................
#...............................................
library("raster")
library("sf")
library("ggplot2")

#...............................................
#...............................................
# define the extent 
e <- extent(-86, -85.3, 12.8, 13.42)

#...............................................
#...............................................
# this file is an illustration for the idea of centralized approach,
# recommendations based on administrative regions
l <- getData("GADM", country = "NIC", level=1)

# crop file for the extent
l <- crop(l, e)

l <- st_as_sf(l)

l$NAME_1 <- c("Var A", "Var A", "Var B")

colors <- c("#ef8a62", "#2166ac","#fee090","#67a9cf")
p <- ggplot(l) +
  geom_sf(aes(fill = as.factor(NAME_1)),
          lwd = 0) + 
  scale_fill_manual(values = colors,
                    labels = c("Var A", "Var B")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 20, colour="black"))

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

d <- crop(d, e)

plot(d)


# replace some values
extf <- extent(c(-86,-85.65,12.8,13.14))
extf2 <- extent(c(-85.8,-85.42,13.25,13.42))
extf3 <- extent(c(-85.5,-85.40,13.13,13.25))
extf4 <- extent(c(-85.7,-85.41,12.8,12.9))
plot(d)
plot(extf, add = TRUE)
plot(extf2, add = TRUE)
plot(extf3, add = TRUE)
plot(extf4, add = TRUE)
# Change values
d[extf][d[extf] == 3] <- 4
d[extf][d[extf] == 1] <- 4

d[extf2][d[extf2] == 2] <- 1

d[extf3][d[extf3] == 1] <- 2
d[extf3][d[extf3] == 3] <- 2

d[extf4][d[extf4] == 1] <- 4
d[extf4][d[extf4] == 3] <- 4

plot(d)

d <- as.data.frame(d, xy = TRUE)

unique(d[,3])

summary(as.factor(d[,3]))

d <- d[!is.na(d$fig1), ]

head(d)

colors <- c("#ef8a62", "#2166ac","#fee090","#67a9cf")

p2 <- 
  ggplot(d) + 
  geom_raster(aes(y = y, x = x, 
                  fill = factor(fig1, levels = c("1","4","3","2")))) +
  scale_fill_manual(values= colors,
                    name = NULL,
                    labels = c("Var A", "Var B", "Var C", "Var D")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 20, colour="black"))

p2
ggsave("output/fig1/decentralized.svg", 
       plot = p2,
       width = 11,
       height = 11, 
       units = "cm")

# remove the layer file
# this can be downloaded again from GADM
file.remove("gadm36_NIC_1_sp.rds")
