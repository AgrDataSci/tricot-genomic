# Model performance and across season forecast
library("tidyverse")
library("magrittr")
library("janitor")
library("svglite")
library("patchwork")
library("gosset")
library("raster")
library("sf")

source("script/helper_00_functions.R")

# ......................................
# .....................................
# Output dir
output <- "output/spatial_dist_kendall_farm_oa/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

# Read data ####
df <- read_csv("data/durumwheat.csv")

load("processing/plmodels/models.rda")

# climatology
load("processing/climatology/climatology.rda")

# environmental indices
ind <- read_csv("data/environmental_indices.csv")

head(ind)

predict(gen$raw$models[[1]])


G <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id")

G <- G[1:length(G),, as.rankings = FALSE]

folds <- gen$raw$folds

kendall <- rep(NA, length(folds))

dt <- gen$raw$data

for(i in seq_along(unique(folds))){
  y <- predict(gen$raw$models[[i]], newdata = dt[folds == i, ], type = "rank")
  x <- G[folds == i, ]
  
  kt <- rep(NA, dim(x)[[1]])
  for(j in seq_len(dim(x)[[1]])){
    
    kt[j] <- kendallTau(y[j, ], x[j, ])[[1]]
    
  }
  print(mean(kt))
  kendall[folds == i] <- mean(kt)
}

r <- data.frame(lon = ind$lon, lat = ind$lat, kendall)

s <- st_as_sf(r, coords = c("lon", "lat"), crs = 4326)
s

eth <- getData("GADM", country = "ETH", level=1)
e <- extent(eth)
e[1] <- 36
e[2] <- 42
e[3] <- 8.5 
e[4] <- 14


eth <- crop(eth, e)

eth <- st_as_sf(eth)

colpall <- colorRampPalette(c("#FFFFFF", "#FFFF80", "#38E009","#1A93AB", "#0C1078"))


ggplot() +
  geom_sf(data = eth, color = alpha("black",0.2)) +
  geom_point(data = r, aes(x = lon, y = lat, color = kendall), size = 0.6, pch = 19) +
  #geom_sf(data = s, aes(color = kendall), size = 0.6, pch = 20) +
  geom_text(aes(label = "Addis Ababa", x = 38.7726, y = 9.0201), size = 3) + 
  theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank()) +
  labs(x="", y="")

ggsave("output/SI/fig_s11_kendall_tau.png",
       plot = last_plot(),
       width = 20,
       height = 20,
       units = "cm")


