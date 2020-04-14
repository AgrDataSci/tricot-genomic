# ..........................................
# ..........................................
# Read and clean durum wheat in Ethiopia
# ..........................................
# ..........................................
library("tidyverse")
library("magrittr")
library("janitor")
#.....................................
#.....................................

# read csv file
mydata <- read_csv("data/raw/climmob_ethiopia_2017.csv", 
                   na = c("#VALUE!", "NA", "<NA>", ""),
                   col_types = cols(GY_gm = col_number()))


mydata %<>% 
  as_tibble(.name_repair = make_clean_names) %>%
  rename(lon = longtiude,
         lat = latitude,
         id  = farmer_no) %>% 
  mutate_if(is.factor, as.character)

# drop columns
drop <- !grepl("code_", names(mydata))

mydata <- mydata[, drop]

drop <- c("no", "altitude", "sex", "days2maturity", "biomass_gm", "genotype")

drop <- !names(mydata) %in% drop

mydata <-  mydata[, drop]

names(mydata)

#.....................................
#.....................................

# change names of accessions
# rename variety Hetosa and Assasa
mydata %<>% 
  mutate(accession = ifelse(accession == "Hetosa",
                            "Hitosa", 
                            accession),
         accession = ifelse(accession == "Assassa",
                            "Asassa", 
                            accession))

#.....................................
#.....................................

# drop all missing rankings or rankings with missing Accession (item id)
mydata %<>%
  filter(!is.na(farmer_rank) & !is.na(accession))

# only keep strict rankings of four distinct items
mydata %>%
  group_by(id) %>%
  summarise(keep = setequal(farmer_rank, 1:4) & n_distinct(accession) == 4) %>%
  filter(keep) ->
  keep

keep <- mydata$id %in% keep$id


mydata %<>% 
  filter(keep)



# check the number of observations per variety
# remove those tested in less than 10 farms 
rmitem <- 
  mydata %>% 
  group_by(accession) %>%  
  count(accession) %>%
  filter(n < 10) %>%
  dplyr::select(accession) %>%
  as.matrix() %>%
  as.vector()


keep <- !mydata$accession %in% rmitem


mydata <- mydata[keep, ]



#.....................................
#.....................................
# remove rankings with less than 2 items
mydata %>%
  group_by(id) %>%
  summarise(keep = n_distinct(accession) >= 2) %>%
  filter(keep) ->
  keep

keep <- mydata$id %in% keep$id


mydata %<>% 
  filter(keep)

#.....................................
#.....................................

# add genotype codes
g <- "data/raw/whoiswho.diversity.panel.txt"
g %<>%  
  read_table2() %>% 
  rename(accession = ID,
         genotype = DNA_CODE) %>% 
  select(accession, genotype) %>% 
  filter(!grepl("_B", genotype))


mydata %<>%
  mutate(accession = tolower(accession))

mydata %<>% 
  merge(. , g, all.x = TRUE, by = "accession") %>% 
  as_tibble() 


mydata %<>% 
  mutate(genotype = ifelse(is.na(genotype), accession, genotype))




#.....................................
#.....................................
# planting dates 

pdate <- read_csv("data/raw/ethiopia_planting_dates.csv")
pdate %<>% 
  as_tibble(.name_repair = janitor::make_clean_names) %>% 
  mutate_if(is.factor, as.character)

mydata <- merge(mydata, pdate[,c("farmer","planting_date")],
                by = "farmer", all.x = TRUE)

mydata$planting_date <- as.Date(mydata$planting_date, "%d/%m/%Y")

mydata$planting_date <- as.integer(mydata$planting_date)

# fill missing planting dates with average per year
sum(is.na(mydata$planting_date))

for(i in seq_along(unique(mydata$year))){
  
  y_i <- unique(mydata$year)[i]
  
  mydata$planting_date <- ifelse(is.na(mydata$planting_date) & mydata$year == y_i, 
                                 mean(mydata$planting_date[mydata$year == y_i], na.rm=TRUE), 
                                 mydata$planting_date)
}

mydata$planting_date <- as.Date(mydata$planting_date, origin = "1970-01-01")

mydata$year <- as.integer(strftime(mydata$planting_date, "%Y"))

sum(is.na(mydata$planting_date))

#.....................................
#.....................................
# lon lat 
sum(is.na(mydata[,c("lon","lat")]))

# debug lat and lon using mean per village 
summary(as.factor(mydata$kebele))

for (j in unique(mydata[,"kebele"])){
  
  mydata[,"lon"] <- ifelse(mydata[,"kebele"] == j & is.na(mydata[,"lon"]),
                           mean(mydata[,"lon"][mydata[,"kebele"] == j ], na.rm = TRUE),
                           mydata[,"lon"])
  
  mydata[,"lat"] <- ifelse(mydata[,"kebele"] == j & is.na(mydata[,"lat"]),
                           mean(mydata[,"lat"][mydata[,"kebele"] == j ], na.rm = TRUE),
                           mydata[,"lat"])
}

sum(is.na(mydata[,c("lon","lat")]))

# debug lat and lon using mean per district
summary(mydata$district)
for (j in unique(mydata[,"district"])){
  mydata[,"lon"] <- ifelse(mydata[,"district"] == j & is.na(mydata[,"lon"]),
                           mean(mydata[,"lon"][mydata[,"district"] == j ], na.rm = TRUE),
                           mydata[,"lon"])
  
  mydata[,"lat"] <- ifelse(mydata[,"district"] == j & is.na(mydata[,"lat"]),
                           mean(mydata[,"lat"][mydata[,"district"] == j ], na.rm = TRUE),
                           mydata[,"lat"])
}
sum(is.na(mydata[,c("lon","lat")]))


# add season 
mydata$season <- paste("Meher", mydata$year, sep="-")
mydata$season <- gsub("20","",mydata$season)


mydata <- as_tibble(mydata)

#.....................................
#.....................................
# check agronomic data

# normalise by plot size
mydata %<>% 
  mutate(gy_gm_raw = gy_gm,
         plot_size = ifelse(year == 2013 & region == "Tigray", 1.6,
                            ifelse(year == 2013 & region != "Tigray", 0.4,
                                   ifelse(year == 2014, 1.6,
                                          ifelse(year == 2015, 1.2, NA)))),
         gy_gm = gy_gm / plot_size)

# mean grain yield
boxplot(mydata$gy_gm ~ mydata$accession, las = 2)

# remove outliers in these variables
y <- split(mydata, mydata$accession)

y <- lapply(y, function(x) {
  i <- c("gy_gm")
  x[,i] <- lapply(x[,i], function(z) {
    
    out <- boxplot.stats(z)$out
    
    z <- ifelse(z %in% out, NA, z)
    
    z
    
  })
  
  x
  
})



y <- bind_rows(y)

boxplot(y$gy_gm ~ y$accession, las = 2)

mydata <- y

boxplot(mydata$gy_gm ~ mydata$accession, las = 2)

names(mydata)
# organize colunm order
mydata %<>% 
  rename(plot_id = plot_no)


first <- c("id","genotype","accession","plot_id","farmer")

mydata <- mydata[c(match(first, names(mydata)),
                   which(!names(mydata) %in% first))] 

mydata %<>% 
  arrange(id)

drop <- c("farmer")

mydata <- mydata[, !names(mydata) %in% drop]

write_csv(mydata, "data/durumwheat.csv")



rev(sort(table(mydata$id)))


