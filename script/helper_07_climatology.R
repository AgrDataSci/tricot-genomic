# Capture trends in climatology over the last 15 years
# across the research area

library("tidyverse")
library("magrittr")
library("janitor")
library("dismo")
library("rgeos")
library("raster")
library("ClimMobTools")
library("gosset")

#.......................................
#.......................................
# Read data ####
# tricot data
df <- read_csv("data/durumwheat.csv")

# station data
load("data/diversity.panel.data.gp.Rdata")
rm(snp.pos, info, geno, farm)

# item names
items <- unique(df$genotype)

df %<>% 
  dplyr::select(id, genotype, lat, lon, planting_date, year)




#......................................
#......................................
# Define time span for phenological stages ####
met %<>% 
  rename(genotype = ID) %>% 
  as_tibble(.name_repair = janitor::make_clean_names) 

keep <- met$genotype %in% items

met <- met[keep, ]

# take the means for each genotype
ts <- 
  met %>% 
  summarise(db = as.integer(median(db, na.rm = TRUE)),
            df = as.integer(median(df, na.rm = TRUE)),
            dm = as.integer(median(dm, na.rm = TRUE)))


span <- 
  ts %>% 
  mutate(veg = db,
         rep = (df - db) + 15,
         gra = dm - df,
         sow2rep = df,
         sow2gra = dm) %>% 
  dplyr::select(-db, -df, -dm)



# # Dates for booting, flowering and maturity
# dates <- 
#   df %>% 
#   select(planting_date, db, df, dm) %>% 
#   mutate(veg = planting_date,
#          rep = (planting_date + db) -8,
#          gra = planting_date + df,
#          sow2rep = planting_date,
#          sow2gra = planting_date) %>% 
#   select(-planting_date, -db, -df, -dm)


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


#......................................
#......................................
# Get planting dates ####

pdates <- 
  df %>% 
  distinct(id, .keep_all = TRUE) %>% 
  dplyr::select(planting_date)
  
pdates %<>% 
  mutate(week = as.integer(format(planting_date, "%U")),
         day = as.integer(format(planting_date,"%j")))

quantile(pdates$week)

years <- 2001:2015
weeks <- 29:31

dates <- integer()

for(i in seq_along(years)) {
  for(j in seq_along(weeks)) {
    d <- as.Date(paste(years[[i]], weeks[[j]], 1, sep="-"), "%Y-%U-%u")
    dates <- c(dates, d)
  }
}

#......................................
#......................................
# Get indices ####

climatology <- list()

for(i in seq_along(dates)) {
  print(i)
  do <- rep(dates[[i]], nrow(coord))
  do <- as.Date(do, origin = "1970-01-01")
  
  h <- temperature(coord,
                   day.one = do,
                   span = ts[[1]])
  
  names(h) <- paste0(names(h), "_veg")
  
  climatology[[i]] <- h
}


output <- "processing/climatology/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

save(climatology, coord, file = paste0(output, "climatology.rda"))

climatology[[1]]

# this is just data visualization
x <- NULL
for(i in seq_along(dates)) {
  y <- data.frame(minDT = climatology[[i]][["minDT_veg"]], 
                  minNT = climatology[[i]][["minNT_veg"]],
                  samp = paste0("S",i))
  x <- rbind(x, y)
}


boxplot(x$minDT ~ x$samp)

boxplot(x$minNT ~ x$samp)



