####################################################
###### Read and clean durum wheat in Ethiopia
# Updated 11Jun2019
####################################################

library("tidyverse")
library("magrittr")
library("janitor")
library("sf")

# ..........................................
# ..........................................
# file with some pre-debuged variables ####
df <- "data/raw/DataSheet_CS_All_FINAL.csv"
df %<>% 
  read.csv(.,
           na.strings = c("NA","", " ", "missing","?","#DIV/0!","#REF!")) %>% 
  as_tibble(.name_repair = janitor::make_clean_names) %>%
  rename(lon = longtiude,
         lat = latitude) %>% 
  mutate_if(is.factor, as.character)


# drop all missing rankings or rankings with missing accession
df %<>%
  filter(!is.na(farmer_rank) & !is.na(accession))

# drop some variables
names(df)

drop <- c("no","altitude","sex","code_accession","days2maturity")

drop <- !names(df) %in% drop

df <- df[drop]

# rename variety Hetosa and Assasa
df %<>% 
  mutate(accession = ifelse(accession == "Hetosa",
                            "Hitosa", 
                            accession),
         accession = ifelse(accession == "Assassa",
                            "Asassa", 
                            accession))
summary(as.factor(df$accession))

# ..........................................
# ..........................................
# Debug lon lat data ####
#lon lat 
sum(is.na(df[,c("lon","lat")]))

df %<>% 
  mutate(lon = ifelse(is.na(lat), NA, lon),
         lat = ifelse(is.na(lon), NA, lat))

# Download country border of ethiopia
eth <-
raster::getData("GADM", country = "ETH", path = "data", level = 0) %>% 
  st_as_sf()

# plot(eth["GID_0"], col = "lightgrey")
# 
# df %>% 
#   select(lon, lat) %>%  
#   filter(!is.na(lon)) %>% 
#   as.matrix() %>% 
#   st_multipoint() %>%
#   st_sfc(crs = st_crs(eth)) -> 
#   lonlat
# 
# plot(eth["GID_0"], col = "lightgrey", reset = FALSE)
# plot(lonlat, col = "blue", cex = 1, 
#      bg = "Steelblue1", pch = 21, add = TRUE)


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

#.........................................
#.........................................
# Keep consistent observations #### 

out <- "data/raw/not_in_genotyping.csv"
out %<>% 
  read_csv() %>% 
  select(accession) %>% 
  t()


out <- !df$accession %in% out

# retain the observations with genotype data
df %<>% 
  filter(out)

# add genotype codes
g <- "data/raw/whoiswho.diversity.panel.txt"
g %<>%  
  read_table2() %>% 
  rename(accession = ID,
         genotype = DNA_CODE) %>% 
  select(accession, genotype) %>% 
  filter(!is.na(genotype)) %>% 
  filter(!grepl("_B", genotype))


df %<>%
  mutate(accession = tolower(accession)) %>% 
  merge(. , g, all.x = TRUE, by = "accession") %>% 
  as_tibble() %>% 
  filter(!is.na(genotype))


# create an id for each plot and remove duplicates
df %<>% 
  mutate(plot_id = paste(genotype, farmer, year, sep = "-"),
         id = as.integer(as.factor(farmer))) %>% 
  arrange(id)

# remove all entries that have duplicates 
df <-
  df[!(duplicated(df$plot_id) | duplicated(df$plot_id, fromLast = TRUE)), ]

df %<>% 
  group_by(id) %>% 
  distinct(genotype, .keep_all = TRUE) %>% 
  distinct(farmer_rank, .keep_all = TRUE) %>% 
  as_tibble()


# remove those tested in less than 10 farms 
df %>% 
  group_by(accession) %>%  
  count(accession) %>%
  filter(n < 10) %>%
  select(accession) %>%
  t() %>%
  as.vector() ->
  drop


drop <- !df$accession %in% drop

df %<>% 
  filter(drop)


# only keep strict rankings of at least 2 distinct items
df %>%
  group_by(id) %>%
  summarise(keep = n_distinct(genotype) >= 2) %>%
  filter(keep) ->
  keep

# apply the logical vector
id <- df$id %in% keep$id

# take rankings
n <- nrow(keep)

# keep selected observations
df <- df[id,]

df %<>% 
  group_by(id) %>% 
  mutate(plot_id = as.integer(as.factor(plot_id))) %>% 
  as_tibble()

# .......................................
# .......................................
# Fix planting dates ####
# planting dates 

df %>% 
  select(sowing_date) %>% 
  mutate(sowing_date = as.character(sowing_date)) ->
  date

date %<>% 
  mutate(sowing_date = gsub("[/]|[.]","-", sowing_date)) %>% 
  separate(., sowing_date, into = c("x","y","z"), sep = "-") %>% 
  mutate_all(as.integer)

summary(as.factor(date$x))
summary(as.factor(date$y))
summary(as.factor(date$z))

# planting dates must be between month 7 and 8
date %<>% 
  mutate(y = ifelse(y < 7 , 8, y),
         y = ifelse(y > 8, 8, y))

summary(as.factor(date$x))
summary(as.factor(date$y))
summary(as.factor(date$z))

date %<>% 
  mutate(planting_date = ifelse(x > 2000, paste(x, y, z, sep = "-"),
                                ifelse(z > 2000, paste(z, y, x, sep = "-"), NA)))

df %<>% 
  mutate(planting_date = date$planting_date,
         planting_date = as.Date(planting_date, format = "%Y-%m-%d"), 
         planting_date = as.integer(planting_date))

summary(as.factor(df$planting_date))


# fill missing planting dates with average per kebele and year
df %<>% 
  mutate(varing = ifelse(year == 2013, -365, 
                         ifelse(year == 2015, 365, 0)))


df %>% 
  select(planting_date, year) %>% 
  filter(!is.na(planting_date), year == 2014) %>% 
  summarise(mean =  mean(planting_date),
            min = min(planting_date),
            max = max(planting_date),
            median = median(planting_date)) ->
  varing


as.integer(rnorm(1, mean = varing$mean, sd = 10), origin = "1970-01-01")

df_split <- split(df, df$id)

df_split <- lapply(df_split, function(X) {
  i <- as.integer(rnorm(1, mean = varing$mean, sd = 10), origin = "1970-01-01")
  i <- i + unique(X$varing)
  
  X$planting_date <- ifelse(is.na(X$planting_date), 
                            i, X$planting_date)
  
  X
})


df <- bind_rows(df_split)


df %<>% 
  mutate(planting_date = as.Date(planting_date, origin = "1970-01-01"),
         year = as.integer(strftime(planting_date, "%Y"))) %>% 
  select(-sowing_date, -varing)

sum(is.na(df$planting_date))

plot(df$planting_date)


# .......................................
# .......................................
# Check agronomic data ####

# set as numeric
df %<>% 
  mutate(gy_gm = as.numeric(gy_gm),
         biomass_gm = as.numeric(biomass_gm),
         mean_seed_no_spike = as.numeric(mean_seed_no_spike),
         mean_tn = as.numeric(mean_tn),
         mean_sl = as.numeric(mean_sl))

# normalise by plot size
df %<>% 
  mutate(gy_gm_raw = gy_gm,
         plot_size = ifelse(year == 2013 & region == "Tigray", 1.6,
                            ifelse(year == 2013 & region != "Tigray", 0.4,
                                   ifelse(year == 2014, 1.6,
                                          ifelse(year == 2015, 1.2, NA)))),
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
  i <- c("mean_tn","mean_seed_no_spike","mean_sl","gy_gm")
  x[,i] <- lapply(x[,i], function(z) {

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

names(df)
# organize colunm order
first <- c("id","genotype","accession","plot_id","farmer")

df <- df[c(match(first, names(df)),
           which(!names(df) %in% first))] 

df %<>% 
  arrange(id)

drop <- c("zone","seedwt_500","biomass_gm","reseason_rank")

df <- df[, !names(df) %in% drop]

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

