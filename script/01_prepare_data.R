####################################################
###### Read and clean durum wheat in Ethiopia
# Updated 06Jun2019
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
         lat = latitude) %>% 
  select(-genotype)

# drop all missing rankings or rankings with missing Accession (item id)
df %<>%
  filter(!is.na(farmer_rank) & !is.na(accession))

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


# ..........................................
# ..........................................
# Drop variables not genotyped ####
out <- "data/raw/not_in_genotyping.csv"
out %<>% 
  read_csv() %>% 
  select(accession) %>% 
  t()


out <- !df$accession %in% out

# retain the observations with genotype data
df %<>% 
  filter(out)

# create an id for each observation
df %<>% 
  mutate(id = paste(plot_no, farmer, sep = "-")) 

# remove all entries that have duplicates 
df <-
  df[!(duplicated(df$id) | duplicated(df$id, fromLast = TRUE)), ]

df %<>% 
  group_by(farmer_no) %>% 
  distinct(accession, .keep_all = TRUE) %>% 
  as_tibble()

# only keep strict rankings of at least 2 distinct items
df %>%
  group_by(farmer_no) %>%
  summarise(keep = n_distinct(accession) >= 2) %>%
  filter(keep) ->
  keep

# apply the logical vector
id <- df$farmer_no %in% keep$farmer_no

# take rankings
n <- nrow(keep)

# keep selected observations
df <- df[id,]


# ..........................................
# ..........................................
# Add plating dates ####
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
  # check which value is the closest to the mean
  # and get the centroid
  z <- mean(df$planting_date[df$year == y_i], na.rm = TRUE)
  z <- unlist(df[df$year == y_i, "planting_date"] - z)
  z <- df[[which.min(abs(z)), "planting_date"]]
  
  df$planting_date <- ifelse(is.na(df$planting_date) & df$year == y_i,
                             z,
                             df$planting_date)

}

df %<>% 
  mutate(planting_date = as.Date(planting_date, origin = "1970-01-01"),
         year = as.integer(strftime(planting_date, "%Y")))


# .......................................
# .......................................
# Debug lon lat data ####
#lon lat 
sum(is.na(df[,c("lon","lat")]))

df %<>% 
  mutate(lon = ifelse(is.na(lat), NA, lon),
         lat = ifelse(is.na(lon), NA, lat))

sum(is.na(df[,c("lon","lat")]))

# Download an outline of ethiopia
raster::getData("GADM", country = "ETH", path = "data", level = 0) %>% 
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


# .......................................
# .......................................
# Dataset with extra variables ####
extra <- "data/raw/DataSheet_CS_All_FINAL.xlsx"
extra %<>% 
  read_excel(., 
             sheet = "data",
             range = cell_cols("A:AA"),
             na = c("NA", " ", "missing","?")) %>% 
  as_tibble(.name_repair = janitor::make_clean_names)


# add id
extra %<>% 
  mutate(id = paste(plot_no, farmer, sep = "-"))

# remove all values with duplicates
# remove all entries that have duplicates 
extra <-
  extra[!(duplicated(extra$id) | duplicated(extra$id, fromLast = TRUE)), ]

# select variables
extra %<>% 
  select(booting_date:id) %>% 
  select(-seedwt_500)

# combine with main dataframe
df %<>% 
  merge(. , extra, by = "id", all.x = TRUE) %>% 
  as_tibble()


# ..........................................
# ..........................................
# Check agronomic data ####

# add plot size
df %<>% 
  mutate(plot_size = ifelse(year == 2013, 0.4, 1.6),
         gy_gm = gy_gm / plot_size)

# mean tiller number
boxplot(df$mean_tn ~ df$accession, las = 2)

# mean seed per spike
boxplot(df$mean_seed_no_spike ~ df$accession, las = 2)

# mean spike length
boxplot(df$mean_sl ~ df$accession, las = 2)

# mean grain yield
boxplot(df$gy_gm ~ df$accession, las = 2)

# remove outliers in these variables
y <- split(df, df$accession)

y <- lapply(y, function(x) {
  
  x[,c(18:21)] <- lapply(x[,c(18:21)], function(z) {

    out <- boxplot.stats(z)$out
    
    z <- ifelse(z %in% out, NA, z)
    
    z
    
  })
  
  x
  
})



y <- bind_rows(y)

boxplot(y$mean_tn ~ y$accession, las = 2)
boxplot(y$mean_seed_no_spike ~ y$accession, las = 2)
boxplot(y$mean_sl ~ y$accession, las = 2)
boxplot(y$gy_gm ~ y$accession, las = 2)

df <- y

boxplot(df$gy_gm ~ df$accession, las = 2)

# id into integers
df %<>% 
  mutate(plot_id = as.integer(as.factor(id)),
         id = as.integer(as.factor(farmer)))


# add genotype codes
gnt <- "data/raw/whoiswho.diversity.panel.txt"
gnt %<>%  
  read_table2() %>% 
  rename(accession = ID,
         genotype = DNA_CODE) %>% 
  select(accession, genotype)

df %<>% 
  merge(. , gnt, all.x = TRUE, by = "accession") %>% 
  as_tibble()


names(df)

# organize colunm order
df <- df[c(2,24,1,23,22,3:21)]

write_csv(df, "data/durumwheat.csv")


# # fix issues in dates ##
# extra %>% 
#   select(id, booting_date, flowering_date, maturity_date) -> 
#   mdate
# 
# mdate[2:4] <- lapply(mdate[2:4], function(x) {
#   
#   check <- gsub("/", "-", x, NA)
#   
#   check <- strsplit(check, split =  "-")
#   
#   check <- do.call("rbind", check)
#   
#   d1 <- paste(check[,3], check[,2], check[,1], sep = "-")
#   d1 <- as.Date(d1, format = "%Y-%m-%d")
#   
#   d2 <- paste(check[,3], check[,1], check[,2], sep = "-")
#   d2 <- as.Date(d2, format = "%Y-%m-%d")
#   
#   check <- as.Date(check, format = "%d/%m/%Y")
#   
#   x <- ifelse(grepl("-", x), x, NA)
#   x <- as.Date(x, format = "%d-%m-%Y")
#   
#   x <-   ifelse(is.na(x),
#                 check, 
#                 x)
#   x <- as.Date(x, origin = "1970-01-01")
#   return(x)
#   
# })
# 
# 
# # calculate days to flowering and to maturity
# df %<>% 
#   mutate(days2flower = flowering_date - planting_date,
#          days2maturity = maturity_date - planting_date)
# 
# plot(df$days2flower)
# plot(df$days2maturity)
# # remove outliers
# df %<>% 
#   mutate(days2flower = ifelse(days2flower < 40 | days2flower > 110, NA, days2flower))
# 
# 
# # try to impute some missing data
# # sum(is.na(df$days2flower))
# # sum(is.na(df$days2maturity))
# # 
# # boxplot(df$days2flower ~ as.factor(df$region))
# df %<>% 
#   mutate(days2flower = ifelse(is.na(days2flower), 
#                               median(days2flower, na.rm = TRUE), 
#                               days2flower),
#          days2maturity = ifelse(is.na( days2maturity), 
#                                 median( days2maturity, na.rm = TRUE), 
#                                 days2maturity))
# plot(df$days2flower)
# plot(df$days2maturity)


# # reshape values
# df %>%
#   select(farmer_no, plot_no, accession) %>%
#   spread(.,
#          key = plot_no,
#          value = "accession") ->
#   vars
# 
# names(vars)[2:5] <- paste0("variety_", letters[1:4])
# 
# df %>%
#   select(farmer_no, plot_no,farmer_rank) %>%
#   spread(.,
#          key = plot_no,
#          value = "farmer_rank") ->
#   R
# 
# names(R)[2:5] <- paste0("rank_variety_", letters[1:4])
# 
# # merge data
# df %<>%
#   select(-accession, -farmer_rank, -plot_no) %>%
#   distinct(farmer_no, .keep_all = TRUE) %>%
#   inner_join(. , vars, by = "farmer_no",all.x = TRUE) %>%
#   inner_join(. , R, by = "farmer_no", all.x = TRUE)

