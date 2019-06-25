# 
library("tidyverse")
library("magrittr")
library("PlackettLuce")
library("gosset")


#...................................
#...................................
# read data 
list.files("data")

# farmer rankings
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()



# additive matrix from genes
load("data/additive.matrix.cs.rda")


# normal prior object from additive matrix
prior <- list(mu = rep(0, ncol(additivemat)),
              Sigma = additivemat)


#.....................................
#.....................................
# farmer rank ####

# create PlackettLuce rankings 
R <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)


mod <- PlackettLuce(R, npseudo = 0)

modprior <- PlackettLuce(R, normal = prior, gamma = TRUE)


AIC(mod)
AIC(modprior)

mod$adherence

modprior$adherence

modprior$normal

