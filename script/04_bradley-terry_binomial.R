library("tidyverse")
library("magrittr")
library("gosset")
library("BradleyTerry2")



df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

output <- "output/exploring/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

#.....................................
#.....................................
# farmer rank ####

# create grouped rankings 
G <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)


object <- to_paircomp(G)




