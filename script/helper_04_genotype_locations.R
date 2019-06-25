



library("tidyverse")
library("magrittr")
library("janitor")
library("sf")

pass <- "data/passport_data_durumwheat.csv"
pass %<>% 
  read_csv() 


eth <-
  raster::getData("GADM", country = "ETH", path = "data", level = 1) %>% 
  st_as_sf()

pass %>%
  dplyr::select(lon, lat) %>%
  filter(!is.na(lon)) ->
  xy
  
xy %>% 
  as.matrix() %>%
  st_multipoint() %>%
  st_sfc(crs = st_crs(eth)) ->
  lonlat

plot(eth["GID_1"], col = "#F2F2F2", reset = FALSE)
plot(lonlat, col = "blue", cex = 1,
     bg = "Steelblue1", pch = 21, add = TRUE)


