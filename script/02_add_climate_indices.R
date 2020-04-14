# Add environmental variables using planting dates 
# and lon lat
library("tidyverse")
library("magrittr")
library("janitor")
library("gosset")
library("climatrends")
library("caret")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/02_add_climate_indices_session_info.txt")


#.......................................
#.......................................
# Read data ####
# tricot data
df <- read_csv("data/durumwheat.csv")

# station data
load("data/diversity.panel.data.gp.rda")
rm(snp.pos, info, geno, farm)

items <- unique(df$genotype)

df %<>% 
  select(id, genotype, lat, lon, planting_date, year)


#......................................
#......................................
# Define time span for phenological stages ####
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


# impute values for missing genotypes
df %<>%
  mutate(db = ifelse(is.na(db), mean(db, na.rm = TRUE), db),
         df = ifelse(is.na(df), mean(df, na.rm = TRUE), df),
         dm = ifelse(is.na(dm), mean(dm, na.rm = TRUE), dm))

ts <-
  df %>% 
  group_by(id) %>% 
  summarise(db = as.integer(max(db, na.rm = TRUE)),
            df = as.integer(max(df, na.rm = TRUE)),
            dm = as.integer(max(dm, na.rm = TRUE)))


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
dates <- 
  df %>% 
  select(planting_date, db, df, dm) %>% 
  mutate(veg = planting_date,
         rep = (planting_date + db) -8,
         gra = planting_date + df,
         sow2rep = planting_date,
         sow2gra = planting_date) %>% 
  select(-planting_date, -db, -df, -dm)

span <- 
  df %>% 
  select(planting_date, db, df, dm) %>% 
  mutate(veg = db,
         rep = (df - db) + 15,
         gra = dm - df,
         sow2rep = df,
         sow2gra = dm) %>% 
  select(-planting_date, -db, -df, -dm)


# ............................................
# ............................................
# # rainfall indices ####
rain <- NULL
for(i in seq_len(ncol(span))) {

  r <- rainfall(df[,c("lon","lat")],
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
  h <- temperature(df[c("lon","lat")],
                   day.one = dates[[i]],
                   span = span[[i]])
  
  names(h) <- paste(names(h), names(dates)[i], sep = "_")
  
  temp <- bind_cols(temp, h)
}


temp

# check for -Inf values
drop <- !apply(temp[1:ncol(temp)], 2, function(x) any(is.infinite(x)))
drop <- as.vector(drop)

temp <- temp[drop]

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
         xy = df$lon + df$lat,
         yx = df$lon - df$lat,
         id = df$id)


write_csv(indices, "data/environmental_indices.csv")

indices
