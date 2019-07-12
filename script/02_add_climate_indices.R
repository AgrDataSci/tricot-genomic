# Add environmental variables using planting dates 
# and lon lat

library("tidyverse")
library("magrittr")
library("janitor")
library("ClimMobTools")
library("caret")

#.......................................
#.......................................
# Read data ####
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

items <- unique(df$genotype)

df %<>% 
  select(id, genotype, lat, lon, planting_date, year)


load("data/chirps.rda")

load("data/modis.rda")

#......................................
#......................................
# Define time span for phenological phases ####
# using station data
load("data/diversity.panel.data.gp.Rdata")

met %<>% 
  rename(genotype = ID) %>% 
  as_tibble(.name_repair = janitor::make_clean_names) 

keep <- met$genotype %in% items

met <- met[keep, ]

# take the means for each genotype
met %<>% 
  group_by(genotype) %>% 
  summarise(db = as.integer(mean(db, na.rm = TRUE)),
            df = as.integer(mean(df, na.rm = TRUE)),
            dm = as.integer(mean(dm, na.rm = TRUE)))


# add it to the main dataset
df %<>% 
  merge(. , met, by = "genotype", all.x = TRUE) %>% 
  as_tibble() %>% 
  arrange(id)


df %>% 
  group_by(id) %>% 
  summarise(db = as.integer(max(db, na.rm = TRUE)),
            df = as.integer(max(df, na.rm = TRUE)),
            dm = as.integer(max(dm, na.rm = TRUE))) ->
  ts


# keep unique id values in the main dataset
df %<>%
  distinct(id, .keep_all = TRUE)

# and combine it with ts
# add it to the main dataset
df %<>% 
  select(-db, -df, -dm) %>% 
  merge(. , ts, by = "id", all.x = TRUE) %>% 
  as_tibble() %>% 
  arrange(id)


# Dates for booting, flowering and maturity
df %>% 
  select(planting_date, db, df, dm) %>% 
  mutate(veg = planting_date,
         rep = (planting_date + db) -8,
         gra = planting_date + df,
         sow2rep = planting_date,
         sow2gra = planting_date) %>% 
  select(-planting_date, -db, -df, -dm) -> 
  dates

df %>% 
  select(planting_date, db, df, dm) %>% 
  mutate(veg = db,
         rep = (df - db) + 15,
         gra = dm - df,
         sow2rep = df,
         sow2gra = dm) %>% 
  select(-planting_date, -db, -df, -dm) -> 
  span


# ............................................
# ............................................
# rainfall indices ####
rain <- NULL
for(i in seq_len(ncol(span))) {
  
  r <- rainfall(chirps,
                day.one = dates[[i]],
                span = span[[i]])
  
  names(r) <- paste(names(r), names(dates)[i], sep = "_")
  
  rain <- bind_cols(rain, r)
}


rain

# check for -Inf values
drop <- !apply(rain[1:ncol(rain)], 2, function(x) any(is.infinite(x)))
drop <- as.vector(drop)

rain <- rain[drop]

# ............................................
# ............................................
# temperature indices ####
temp <- NULL
for(i in seq_len(ncol(span))) {
  h <- temperature(modis_approx,
                   day.one = dates[[i]],
                   span = span[[i]])
  
  names(h) <- paste(names(h), names(dates)[i], sep = "_")
  
  temp <- bind_cols(temp, h)
}


temp


# combine temperature and rainfall indices
indices <- bind_cols(temp, rain)

# remove variables with low variance
out <- nearZeroVar(indices, freqCut = 95/50, uniqueCut = 50)

print(names(indices)[out])

indices <- indices[-out]

indices %<>% 
  mutate(year = df$year,
         lat = df$lat,
         lon = df$lon,
         id = df$id)


write_csv(indices, "data/environmental_indices.csv")


indices
