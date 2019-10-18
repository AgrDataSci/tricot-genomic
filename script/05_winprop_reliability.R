library("tidyverse")
library("magrittr")
library("dismo")
library("rgeos")
library("raster")
library("gosset")



#.......................................
#.......................................
# Read data ####
# tricot data
df <- read_csv("data/durumwheat.csv")

# climatology
load("processing/climatology/climatology.rda")

# PL models
load("processing/plmodels/models.rda")

# yield from station data
load("data/diversity.panel.data.gp.Rdata")
rm(snp.pos, geno, info, farm)

# genotypes' names
items <- unique(df$genotype)[unique(df$genotype) %in% unique(met$ID)]

# names of improved varieties recommended by the Wheat Atlas
improved <- c("asassa", "ude", "hitosa")
# genotypes of improved varieties
improved <- unique(df[df$accession %in% improved, "genotype"])[[1]]

#......................................
#......................................
# Sampled points ####

# coordinates data
coord <- df[,c("lon","lat")]

# create a convexHull to limit the model to the 
# areas where trials where evaluated
# calculate largest distance
largedist <- 
  coord %>%
  raster::pointDistance(., longlat = FALSE) %>%
  max(., na.rm = TRUE)

# make a convex hull 
hull <- convHull(coord, lonlat = TRUE)

# extent convex hull
hull <- gBuffer(hull@polygons, 
                width = 0.15 * largedist)

crs(hull) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# convert into a raster
r <- raster(hull, res = 0.0416, ext = extent(hull))

hull <- rasterize(hull, r, field = 1,
                  background = NA)


coord <- as.data.frame(hull, xy = TRUE)

coord <- coord[!is.na(coord[,3]), -3 ]

set.seed(345)
s <- sample(1:nrow(coord), 1200)

coord <- coord[s, ]


# ................................
# ................................
# predictions

names(coord) <- c("lon","lat")
coord %<>% 
  mutate(xy = lon + lat,
         yx = lon - lat)

loc_pred <- predict(loc, newdata = coord)


env_gy <- array(NA, 
                dim = c(nrow(coord),
                        length(items),
                        length(climatology)),
                dimnames = list(1:nrow(coord),
                                items,
                                1:length(climatology)))
for(i in seq_along(climatology)){
  print(i)
  y <- worst_regret(climatology[[i]]))
  env_gy[,items,i] <- y[,items]
}


env_gy <- apply(env_gy, c(1,2), mean)


object <- env_gy
