# Find matches between citizen science dataset
# and a farmer survey in Ethiopia
#
# Kauê de Sousa

library(tidyverse)
library(magrittr)

# read Martina's household data with her ids
hh <- "data/HH_id2019_Ethiopia.csv"
hh %<>% 
  read_csv()

# also the file with kebeles ids
kb <- "data/HH_id2019_Ethiopia_kebeles.csv"
kb %<>% 
  read_csv() %>% 
  select(-district_id)

# join the two datasets
hh %<>% 
  inner_join(., kb, by = "kebele_id")

# read the citizen science data
df <- "data/durumwheat.csv"
df %<>% 
  read_csv() %>% 
  select(id, farmer, region, kebele) %>% 
  distinct(id, .keep_all = TRUE)

# get a vector with unique names for kebeles
kebele <- sort(unique(df$kebele))

kebele[!kebele %in% sort(unique(hh$kebele))]

# fix some issues in kebeles names
# first colunm is the names that appear in the citizen science data
# the second is in Martina's survey
repl <- matrix(c("Adi kunti","Adi Kuenti",
                 "Koka","Koca",
                 "Burtilik","Birtilik",
                 "Medagolat","Meda golat" ,
                 "Mysedri","May sedri",
                 "Tility","Tlity",
                 "Weketa","Wekete",
                 "Werkaye","Workaye",
                 "Wole","Wole Amba",
                 "Weyn amba","Weyra Amba"), 
               nrow = 10, ncol = 2, byrow = TRUE)

# replace names in the citizen science dataset
for(i in seq_along(repl[,1])) {
  df$kebele <- ifelse(df$kebele == repl[i,1],
                      repl[i,2],
                      df$kebele)
}

# get a new vector with kebele names
kebele <- sort(unique(df$kebele))

# check the kebeles in Martina's data that appear in the 
# citizen science dataset
kebele <- kebele[kebele %in% sort(unique(hh$kebele))]



# keep only kebeles that occurs in both datasets
keep <- df$kebele %in% kebele
df <- df[keep, ]

keep <- hh$kebele %in% kebele
hh <- hh[keep,]

# empty vector to match with survey dataset
df$hhid <- NA
df$hh_member <- NA

# run over kebeles 
# then over farmers names within kebeles to find possible match
for (j in seq_along(kebele)) {
  # subset datasets
  x <- hh[hh$kebele == kebele[j], ]

  y <- df[df$kebele == kebele[j], ]

  for(i in seq_along(y$id)) {
    # find index for possible matches
    # we can play arround these arguments and try
    # to increase number of matches
    index <- agrep(y[[i, "farmer"]], x$hh_member, max = 3, ignore.case = TRUE)
    
    # if no match, go to the next farmer
    if(length(index) == 0) next
    
    # if more than one match, take the first
    if(length(index) > 1) index <- index[1]
    
    # get the hh id
    index_id <- x$hhid[index]
    
    # and the farmer name
    index_name <- x$hh_member[index]
    
    # put it into the citizen science dataset
    df[[df$id == y$id[i], "hhid"]] <- index_id
    
    df[[df$id == y$id[i], "hh_member"]] <- index_name

  }

}

df

sum(!is.na(df$hhid))







