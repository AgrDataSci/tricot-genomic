# Model scenarios with explanatory variables 
# to predict variety performance using PlackettLuce models

library("tidyverse")
library("magrittr")
library("PlackettLuce")
library("gosset")
library("foreach")
library("doParallel")
library("abind")


#.................................................
#.................................................
# Data ####
list.files("data")

# farmer rankings
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()


# additive matrix from gene markers
load("data/additive.matrix.cs.rda")


# environmental indices
ind <- "data/environmental_indices.csv"
ind %<>% 
  read.csv()

#.................................................
#.................................................
# PlackettLuce rankings ####
G <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)

# normal prior coefficients
mod <- PlackettLuce(G)
prior <- list(mu = as.vector(coef(mod)),
              Sigma = additivemat)

items <- names(coef(mod))


#.................................................
#.................................................
# Forward selection ####

output <- "output/plmodels/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

# cross validation parameters
n <- length(G)
minsize <- floor(n*0.2)
gamma <- TRUE
alpha <- 0.01
npseudo <- 5

set.seed(123)
k <- n
folds <- sample(rep(1:k, times = ceiling(n / k), length.out = n))
#folds <- as.integer(as.factor(ind$year))
#k <- max(folds)
keep <- grepl("DT_|NT_|DTR|year|lon|lat", names(ind))

data <- cbind(G, empty_model = 0, ind[keep])

#.................................................
#.................................................
# Model scenarios ####
cvnull <- crossvalidation(G ~ empty_model,
                          data = data,
                          k = k,
                          folds = folds,
                          npseudo = npseudo,
                          alpha = alpha)

cvclim <- crossvalidation(G ~ minNT_sow2gra + minDT_veg,
                          data = data,
                          k = k,
                          folds = folds,
                          gamma = gamma,
                          minsize = minsize,
                          npseudo = npseudo,
                          alpha = alpha)

cvclimgen <- crossvalidation(G ~ minNT_sow2gra + minDT_veg,
                          data = data,
                          k = k,
                          folds = folds,
                          gamma = gamma,
                          normal = prior,
                          minsize = minsize,
                          npseudo = npseudo,
                          alpha = alpha)


cvnull
cvclim
cvclimgen


capture.output(print(cvnull),
                print(cvclim),
                print(cvclimgen),
                file = paste0(output, "model_estimates.txt"))
 
# null <- pltree(G ~ empty_model,
#                data = data,
#                npseudo = npseudo,
#                alpha = alpha)
# 
# design <- pltree(G ~ lon + lat + year,
#                data = data,
#                minsize = minsize,
#                npseudo = npseudo,
#                alpha = alpha)
# 
# clim <- pltree(G ~ minNT_sow2gra + minDT_veg,
#                data = data,
#                gamma = gamma,
#                minsize = minsize,
#                npseudo = npseudo,
#                alpha = alpha)
# 
# climgen <- pltree(G ~ minNT_sow2gra + minDT_veg,
#                   data = data,
#                   gamma = gamma,
#                   normal = prior,
#                   minsize = minsize,
#                   npseudo = npseudo,
#                   alpha = alpha)
# 
# pseudoR2(null)
# pseudoR2(design)
# pseudoR2(clim)
# pseudoR2(climgen)
# 
# str(null)
# 
# 
