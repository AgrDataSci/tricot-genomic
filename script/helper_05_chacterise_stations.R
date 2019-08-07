# characterise climate in station plots 
library("ClimMobTools")

lonlat <- matrix(c(38.866667 , 11.666667,
                   39.166667 , 13.650000), 
                 ncol = 2, nrow = 2, byrow = TRUE)


lonlat <- rbind(lonlat, lonlat)

dimnames(lonlat)[[2]]<-c("lon","lat")

pdates <- rep(c("2012-07-05","2013-07-08"), each = 2)
pdates <- as.Date(pdates, format = "%Y-%m-%d")


#......................................
#......................................
# Define time span for phenological phases ####
# using station data
load("data/diversity.panel.data.gp.Rdata")

met %<>% 
  rename(genotype = ID) %>% 
  as_tibble(.name_repair = janitor::make_clean_names) 

# take the means for each genotype
met %<>% 
  group_by(year, location) %>% 
  summarise(db = as.integer(mean(db, na.rm = TRUE)),
            df = as.integer(mean(df, na.rm = TRUE)),
            dm = as.integer(mean(dm, na.rm = TRUE))) %>% 
  ungroup()


df <- bind_cols(met, planting_date = pdates)

# calculate spans

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
# # rainfall indices ####
rain <- NULL
for(i in seq_len(ncol(span))) {
  
  r <- rainfall(lonlat,
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
  h <- temperature(lonlat,
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

indices <- cbind(indices, lonlat, df[,c("location","year")])

write_csv(indices, "data/environmental_indices_station.csv")


