library("tidyverse")
library("magrittr")
library("PlackettLuce")
library("gosset")
library("foreach")
library("doParallel")
library("abind")

source("script/helper_07_forward.R")

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


# environmental indices
ind <- "data/environmental_indices.csv"
ind %<>% 
  read.csv()

#.....................................
#.....................................
# create PlackettLuce rankings ####
G <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)

mod <- PlackettLuce(G)

# normal prior object from additive matrix
prior <- list(mu = as.vector(coef(mod)),
              Sigma = additivemat)

#.................................................
#.................................................
# Set parameters for forward selection ####
#out <- !grepl("Rx1day|Rx5day|SU|Rtotal|SDII", names(ind))
#ind <- ind[out]
keep <- grepl("DT_|NT_|lon|lat|year", names(ind))

data <- cbind(G, empty_model = 0, ind[keep])

n <- nrow(data)

minsize <- floor(n*0.15)

null <- pltree(G ~ empty_model,
               data = data,
               npseudo = 5,
               normal = prior)

design <- pltree(G ~ lon + lat + year,
               data = data,
               npseudo = 5,
               minsize = minsize)

clim <- pltree(G ~ minDT_veg + minNT_sow2gra,
               data = data,
               npseudo = 5,
               gamma = TRUE,
               minsize = minsize)

climgen <- pltree(G ~ minDT_veg + minNT_sow2gra,
                  data = data,
                  npseudo = 5,
                  gamma = TRUE,
                  normal = prior,
                  minsize = minsize)

predict(null)
predict(design)
predict(clim)
predict(climgen)

pseudoR2(null)
pseudoR2(design)
pseudoR2(clim)
pseudoR2(climgen)



# cross validation 
set.seed(123)
k <- n
folds <- sample(rep(1:k, times = ceiling(n / k), length.out = n))

cvnull <- crossvalidation(G ~ empty_model,
                          data = data,
                          k = k, 
                          folds = folds,
                          npseudo = 5,
                          normal = prior)

# cvclim <- crossvalidation(G ~ minDT_veg + minNT_sow2gra,
#                           data = data,
#                           k = k, 
#                           folds = folds,
#                           gamma = TRUE,
#                           npseudo = 5)

cvclim <- rep(NA, times = k)
for(i in seq_len(k)){
  cat(i)
  train <- data[-folds[i], ]
  test <- data[folds[i], ]
  
  fit <- pltree(G ~ minDT_veg + minNT_sow2gra,
                data = train,
                npseudo = 5,
                gamma = TRUE,
                minsize = minsize)
  
  cvclim[i] <- deviance(fit, newdata = test)
 
}



cvclimgen <- crossvalidation(G ~ minDT_veg + minNT_sow2gra,
                          data = data,
                          k = k, 
                          folds = folds,
                          gamma = TRUE,
                          npseudo = 5,
                          normal = prior)



