# Write maps used in the illustration of Fig 01.
#...............................................
#...............................................

library("raster")
library("sf")
library("ggplot2")
library("svglite")

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

p <- ggplot(l) +
  geom_sf(aes(fill = as.factor(l$NAME_1)),
          lwd = 0) + 
  scale_fill_manual(values = c("#2c7bb6", "#fee090"),
                    labels = c("Var A", "Var B")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 20, colour="black"))


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

d <- as.data.frame(d, xy = TRUE)

d <- d[!is.na(d$fig1), ]

head(d)

p2 <- 
  ggplot(d) + 
  geom_raster(aes(y = y, x = x, 
                  fill = factor(fig1, levels = c("3","2","1")))) +
  scale_fill_manual(values= c("#2c7bb6", "#fee090", "#d73027"),
                    name = NULL,
                    labels = c("Var A", "Var B", "Var C")) +
  theme_void() +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text= element_text(size = 20, colour="black"))


ggsave("output/fig1/decentralized.svg", 
       plot = p2,
       width = 11,
       height = 11, 
       units = "cm")

# remove the layer file
# this can be downloaded again from GADM
file.remove("gadm36_NIC_1_sp.rds")
