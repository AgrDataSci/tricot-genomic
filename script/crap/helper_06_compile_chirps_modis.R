# Add environmental variables using planting dates 
# and lon lat

library("tidyverse")
library("magrittr")
library("ClimMobTools")
library("raster")


#.......................................
#.......................................
# Read data ####
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

df %<>%
  distinct(id, .keep_all = TRUE)

n <- nrow(df)

# define some parameters

# the external drive with MODIS and CHIRPS data
edrive <- "E:/rasters/"

my_proj <- "+proj=longlat +datum=WGS84"

ext <- raster::extent(68, 86, 9, 33)

modis_factor <- 0.02

k_to_c <- 273.15

tiles <- c("h21v07","h21v08")

years <- c(2012:2016)

dates <- NULL

# generate a sequency of 8-days from 2013 to 2015
for(i in c(2012:2015)){
  dates <- c(dates, c(seq(as.Date(paste0(i, "-01-01"), "%Y-%m-%d"), 
                          as.Date(paste0(i, "-12-31"), "%Y-%m-%d"), 8)))
}

dates <- as.Date(c(dates, seq(as.Date("2016-01-01", "%Y-%m-%d"), 
                              as.Date("2016-12-26", "%Y-%m-%d"), 8)), 
                 origin = "1970-01-01")

# and a vector with full dates (1 days)
full_dates <- seq(dates[1], dates[length(dates)], 1)

#.......................................
#.......................................
# Add MODIS data ####

# create an array with 2 layers
# modis day and night
modis <- array(NA, dim = c(n, length(dates), 2), 
               dimnames = list(NULL, as.character(dates), c(1:2)))

# files in this directory has different extensions, 
# so each tile is loaded individually
# they are also very large so it is retrieved from a external drive
# check if the drive exists

for (i in 1:2){
  
  if(i == 1 ) m <- "day" else m <- "night"
  
  for(j in 1:length(tiles)){
    
    file <- list.files(paste0(edrive, "modis/smooth_Ethiopia/", tiles[j]), 
                       pattern = m, 
                       full.names = TRUE)
    
    
    file <- file[grepl(paste0(years, collapse = "|"), file)]
    
    
    file <- stack(file)
    
    x <- extract(file, df[, c("lon", "lat")])
    
    # convert to Celsius
    x <- x * modis_factor - k_to_c
    
    # remove wrong data
    x[x < 1] <- NA
    
    dimnames(x)[[1]] <- 1:n
    
    x <- x[!is.na(x[, 3]),]
    
    modis[as.integer(dimnames(x)[[1]]) , , i] <- x
  
  }
}

rm(x)
sum(is.na(modis))

# replace NA's in temp using the median per day
for (i in 1:dim(modis)[2]) {
  modis[, i, 1] <- ifelse(is.na(modis[, i, 1]),
                          median(modis[, i, 1], na.rm = TRUE),
                          modis[, i, 1])
  modis[, i, 2] <- ifelse(is.na(modis[, i, 2]),
                          median(modis[, i, 2], na.rm = TRUE),
                          modis[, i, 2])
}

# apply linear interpolation to the 8-day temperature
modis_approx <- array(NA, c(n, length(full_dates), 2), 
                      dimnames = list(df$id,  as.character(full_dates), c(1:2) ))

for(i in 1:dim(modis)[1]){
  modis_approx[i,,1] <- approx(modis[i,,1], n = length(full_dates))[[2]]
  modis_approx[i,,2] <- approx(modis[i,,2], n = length(full_dates))[[2]]
}
#save(modis, modis_approx, file = "data/modis.rda")

# #add CHIRPS data
cat("Readind CHIRPS data. Time:", date(), "\n")
#source CHIRPS (precipitation data)
#http://chg.geog.ucsb.edu/data/chirps/
years <- c(2013:2016)
chirps <- list()
#run CHIRPS over the chosen years
for (i in years){
  print(i)
  info <- raster::stack(list.files(paste0(edrive, "CHIRPS/", i), 
                                   pattern = ".tif", full.names = TRUE))
  names(info) <- gsub("chirps.v2.0.", "d", names(info))
  chirps[paste("y", i ,sep="")] <- info
  rm(info)
}
#convert list into one large raster stack
chirps <- raster::stack(chirps)
#chirps <- raster::crop(chirps, my_ext)
#extract precipitation info using lat and lon of each farm
chirps <- raster::extract(chirps, df[,c("lon","lat")])
sum(is.na(chirps))
dimnames(chirps)[[2]] <- gsub("d", "", dimnames(chirps)[[2]])
dimnames(chirps)[[2]] <- gsub("[.]", "-", dimnames(chirps)[[2]])
#save(chirps, file = "data/chirps.rda")



