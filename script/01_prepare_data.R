####################################################
###### Read and clean durum wheat in Ethiopia
# Updated 31May2019
####################################################

library("tidyverse")
library("readxl")
library("magrittr")
library("janitor")
library("sf")

# ..........................................
# ..........................................
# file with some pre-debuged variables ####
df <- "data/raw/climmob_ethiopia_2017.csv"
df %<>% 
  read_csv(., 
           na = c("#VALUE!", "NA", "<NA>", ""))

# tidy up the colnames
df %<>% 
  as_tibble(.name_repair = janitor::make_clean_names) %>%
  rename(lon = longtiude,
         lat = latitude)

# drop all missing rankings or rankings with missing Accession (item id)
df %<>%
  filter(!is.na(farmer_rank) & !is.na(accession))

# only keep strict rankings of at least 2 distinct items
df %>%
  group_by(farmer_no) %>%
  summarise(keep = setequal(farmer_rank, 1:4) & n_distinct(accession) == 4) %>%
  filter(keep) ->
  keep

# apply the logical vector
id <- df$farmer_no %in% keep$farmer_no

# take rankings
n <- nrow(keep)

itemnames <- sort(unique(df$accession[id]))

# keep selected observations
df <- df[id,]

# drop some variables
names(df)
drop <- c("no" , "code_region" , "code_year" , "code_zone" , "code_district",
          "code_kebelle" , "altitude" , "code_farmer" , "sex" , "code_accession",
          "days2maturity" , "gy_gm" , "biomass_gm")

drop <- !names(df) %in% drop

df <- df[drop]


# rename variety Hetosa
df %<>% 
  mutate(accession = ifelse(accession == "Hetosa",
                            "Hitosa", 
                            accession))

# reshape values
df %>% 
  select(farmer_no, plot_no, accession) %>% 
  spread(.,
         key = plot_no,
         value = "accession") ->
  vars 
  
names(vars)[2:5] <- paste0("variety_", letters[1:4])

df %>% 
  select(farmer_no, plot_no,farmer_rank) %>% 
  spread(.,
         key = plot_no,
         value = "farmer_rank") ->
  R 

names(R)[2:5] <- paste0("rank_variety_", letters[1:4])

# merge data
df %<>% 
  select(-accession, -farmer_rank, -plot_no) %>% 
  distinct(farmer_no, .keep_all= TRUE) %>% 
  inner_join(. , vars, by = "farmer_no",all.x = TRUE) %>% 
  inner_join(. , R, by = "farmer_no", all.x = TRUE)


# ..........................................
# ..........................................
# add plating dates ####
# planting dates 
pdate <- "data/raw/ethiopia_planting_dates.csv"
pdate %<>% 
  read_csv() %>% 
  as_tibble(.name_repair = janitor::make_clean_names)

df %<>% 
  merge(., pdate[,c("farmer","planting_date")], 
        by = "farmer", all.x=TRUE) %>% 
  as_tibble()
  
df %<>% 
  mutate(planting_date = as.Date(planting_date, "%d/%m/%Y"),
         planting_date = as.integer(planting_date))

# fill missing planting dates with average per year
for(i in seq_along(unique(df$year))){
  y_i <- unique(df$year)[i]
  df$planting_date <- ifelse(is.na(df$planting_date) & df$year == y_i, 
                                 median(df$planting_date[df$year == y_i], na.rm=TRUE), 
                             df$planting_date)
}

df %<>% 
  mutate(planting_date = as.Date(planting_date, origin = "1970-01-01"),
         year = as.integer(strftime(planting_date, "%Y")))

# .......................................
# .......................................
# add maturity date ####
mdate <- "data/raw/DataSheet_CS_Yoha_All_October_New.xlsx"
mdate %<>% 
  read_excel(., 
             sheet = "data",
             range = cell_cols("A:AA"),
             na = c("NA", " ", "missing")) %>% 
  as_tibble(.name_repair = janitor::make_clean_names)

names(mdate)

mdate %<>% 
  select(plot_no, farmer, flowering_date, maturity_date)

mdate

# create a new variable "check" with the dates that are in the different formart
mdate %<>% 
  mutate(check = ifelse(grepl("[/]", flowering_date), flowering_date, NA),
         flowering_date = ifelse(grepl("-", flowering_date), flowering_date, NA))

# then convert it into dates
mdate %<>% 
  mutate(check = as.Date(check, format = "%d/%m/%Y"),
         flowering_date = as.Date(flowering_date, format = "%d-%m-%Y"))

# and add all into flowering_date
mdate %<>% 
  mutate(flowering_date = ifelse(is.na(flowering_date),
                                 check, flowering_date),
         flowering_date = as.Date(flowering_date, origin = "1970-01-01")) %>% 
  select(-check)


# same for maturity_date
# create a new variable "check" with the dates that are in the different formart
mdate %<>% 
  mutate(check = ifelse(grepl("[/]", maturity_date), maturity_date, NA),
         maturity_date = ifelse(grepl("-", maturity_date), maturity_date, NA))

# then convert it into dates
mdate %<>% 
  mutate(check = as.Date(check, format = "%d/%m/%Y"),
         maturity_date = as.Date(maturity_date, format = "%d-%m-%Y"))

# and add all into maturity_date
mdate %<>% 
  mutate(maturity_date = ifelse(is.na(maturity_date),
                                 check, maturity_date),
         maturity_date = as.Date(maturity_date, origin = "1970-01-01")) %>% 
  select(-check)


# take the median of these dates
mdate %<>% 
  group_by(farmer) %>% 
  summarise(flowering_date = median(flowering_date, na.rm = TRUE),
            maturity_date = median(maturity_date, na.rm = TRUE))

# combine with main dataframe
df %<>% 
  merge(. , mdate, by = "farmer", all.x = TRUE) %>% 
  as_tibble()


# calculate days to flowering and to maturity
df %<>% 
  mutate(days2flower = flowering_date - planting_date,
         days2maturity = maturity_date - planting_date)

plot(df$days2flower)
plot(df$days2maturity)
# remove outliers
df %<>% 
  mutate(days2flower = ifelse(days2flower < 40 | days2flower > 110, NA, days2flower))


# try to impute some missing data
# sum(is.na(df$days2flower))
# sum(is.na(df$days2maturity))
# 
# boxplot(df$days2flower ~ as.factor(df$region))
df %<>% 
  mutate(days2flower = ifelse(is.na(days2flower), 
                              median(days2flower, na.rm = TRUE), 
                              days2flower),
         days2maturity = ifelse(is.na( days2maturity), 
                                median( days2maturity, na.rm = TRUE), 
                                days2maturity))
plot(df$days2flower)
plot(df$days2maturity)

# .......................................
# .......................................
# debug lon lat data ####

#lon lat 
sum(is.na(df[,c("lon","lat")]))

df %<>% 
  mutate(lon = ifelse(is.na(lat), NA, lon),
         lat = ifelse(is.na(lon), NA, lat))

sum(is.na(df[,c("lon","lat")]))

# Let's download an outline of Denmark
raster::getData("GADM", country = "ETH", level = 0) %>% 
  st_as_sf() -> eth

plot(eth["GID_0"], col = "lightgrey")

df %>% 
  select(lon, lat) %>%  
  filter(!is.na(lon)) %>% 
  as.matrix() %>% 
  st_multipoint() %>%
  st_sfc(crs = st_crs(eth)) -> 
  lonlat

plot(eth["GID_0"], col = "lightgrey", reset = FALSE)
plot(lonlat, col = "blue", cex = 1, 
     bg = "Steelblue1", pch = 21, add = TRUE)


# fill NAs using comunity centroid 
summary(as.factor(df$kebele))
varing <- as.vector(unique(df[["kebele"]]))

for (j in seq_along(varing)){
  
  print(varing[j])
  
  df %>% 
    filter(kebele == varing[j]) %>% 
    filter(!is.na(lon)) %>% 
    mutate(z = lon + lat,
           z = z - mean(z)) %>% 
    select(lon, lat, z) ->
    xyz
  
  # check which value is the closest to the mean
  # and get the centroid
  z <- xyz[which.min(abs(xyz$z)), c("lon","lat")]
  
  if(nrow(z) == 0) next
  
  df %<>% 
    mutate(lon = ifelse(kebele == varing[j] & is.na(lon),
                        z[[1,"lon"]],
                        lon),
           lat = ifelse(kebele == varing[j] & is.na(lat),
                        z[[1,"lat"]],
                        lat))
}

sum(is.na(df[,c("lon","lat")]))

# debug lat and lon using mean per district
varing <- as.vector(unique(df[["district"]]))
for (j in seq_along(varing)){
  
  print(varing[j])
  
  df %>% 
    filter(district == varing[j]) %>% 
    filter(!is.na(lon)) %>% 
    mutate(z = lon + lat,
           z = z - mean(z)) %>% 
    select(lon, lat, z) ->
    xyz
  
  # check which value is the closest to the mean
  # and get the centroid
  z <- xyz[which.min(abs(xyz$z)), c("lon","lat")]
  
  if(nrow(z) == 0) next
  
  df %<>% 
    mutate(lon = ifelse(district == varing[j] & is.na(lon),
                        z[[1,"lon"]],
                        lon),
           lat = ifelse(district == varing[j] & is.na(lat),
                        z[[1,"lat"]],
                        lat))
}

sum(is.na(df[,c("lon","lat")]))

df %>% 
  select(lon, lat) %>%  
  filter(!is.na(lon)) %>% 
  as.matrix() %>% 
  st_multipoint() %>%
  st_sfc(crs = st_crs(eth)) -> 
  lonlat

plot(eth["GID_0"], col = "lightgrey", reset = FALSE)
plot(lonlat, col = "blue", cex = 1, 
     bg = "Steelblue1", pch = 21, add = TRUE)



write_csv(df, "data/durumwheat.csv")



