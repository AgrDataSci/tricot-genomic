# Add environmental variables using planting dates 
# and lon lat

library("tidyverse")
library("magrittr")
library("ClimMobTools")
library("caret")

#.......................................
#.......................................
# Read data ####
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

df %<>%
  distinct(id, .keep_all = TRUE)

n <- nrow(df)

load("data/modis.rda")

load("data/chirps.rda")

df


# Average dates for booting, flowering and maturity
# c(60,75,135)
spans <- as.matrix(cbind(period = c("veg", "rep", "gra", "sow2rep", "sow2gra"),
                   ts = c(65, 30, 45, 95, 140),
                   sumtosow = c(0, 65, 95, 0, 0)))


# ............................................
# ............................................
# rainfall indices ####
rain <- NULL
for(i in seq_along(spans[,1])) {
  r <- rainfall(chirps,
                day.one = (df$planting_date + as.integer(spans[i, "sumtosow"])),
                span = as.integer(spans[i, "ts"]))
  
  names(r) <- paste(names(r), spans[i, "period"], sep = "_")
  
  rain <- bind_cols(rain, r)
}


rain




# ............................................
# ............................................
# temperature indices ####
temp <- NULL
for(i in seq_along(spans[,1])) {
  h <- temperature(modis_approx,
                   day.one = (df$planting_date + as.integer(spans[i, "sumtosow"])),
                   span = as.integer(spans[i, "ts"]))
  
  names(h) <- paste(names(h), spans[i, "period"], sep = "_")
  
  temp <- bind_cols(temp, h)
}


temp



# combine temperature and rainfall indices
indices <- bind_cols(rain, temp)

indices %<>% 
  mutate(year = df$year)

# remove variables with low variance
out <- nearZeroVar(indices, freqCut = 95/50, uniqueCut = 50)

names(indices)[out]

indices <- indices[-out]


write_csv(indices, "data/environmental_indices.csv")

