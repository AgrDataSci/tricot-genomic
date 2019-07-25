# create a SpatialPoints object from accessions locations
library("tidyverse")
library("magrittr")
library("rgeos")
library("rgdal")

pass <- "data/passport_data_durumwheat.csv"
pass %<>% 
  read_csv() 

pass %<>% 
  filter(!is.na(lon)) %>% 
  select(genotype, lon, lat)

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

p <- SpatialPointsDataFrame(pass[c("lon","lat")], pass, 
                            proj4string = CRS(proj))



