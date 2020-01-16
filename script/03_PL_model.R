# Model scenarios with explanatory variables 
# to predict variety performance using PlackettLuce models

library("tidyverse")
library("PlackettLuce")
library("gosset")
library("caret")

#.................................................
#.................................................
# Data ####
list.files("data")

# farmer rankings
df <- read_csv("data/durumwheat.csv")

# environmental indices
ind <- read_csv("data/environmental_indices.csv")

# genotypic data
load("data/genotypic.data.durum.wheat.rda")
gene <- imputedcl
rm(imputedcl, hmp)

#.................................................
#.................................................
# PlackettLuce rankings ####
G <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id", 
                  ascending = TRUE,
                  group = TRUE)

#.................................................
#.................................................
# normal prior coefficients ####
mod <- PlackettLuce(G)

items <- names(coef(mod))

genitems <- items[items %in% dimnames(gene)[[1]]]

# keep only the items used in this analysis
gene <- gene[dimnames(gene)[[1]] %in% genitems, ]

# remove those with no variance
drop <- nearZeroVar(gene)
gene <- gene[, -drop]

set.seed(765)
am <- matrix(runif(73*ncol(gene), -1e-1, 1e-1), 
             nrow = length(items), 
             ncol = dim(gene)[[2]],
             dimnames = list(items,
                             dimnames(gene)[[2]]))

am[1:10, 1:10]

am[genitems,] <- gene

am[1:10, ]

am <- t(am)

am <- cov(am)

chol.default(am)

prior <- list(mu = as.vector(coef(mod)),
              Sigma = am)

#.................................................
#.................................................
# Forward selection ####
output <- "processing/plmodels/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

# cross validation parameters
n <- length(G)
minsize <- round(n * 0.30, -2)
npseudo <- 5
alpha <- 0.01
bonferroni <- TRUE
gamma <- TRUE
folds <- as.numeric(as.factor(ind$year))
k <- max(folds)


# select variables
keep <- which(grepl("DT_|NT_|MLDS_rep|lon|lat|xy|yx", names(ind)))

ind <- ind[, keep]

# put all together as mydata
mydata <- cbind(G, empty_model = rep(0, nrow(ind)), ind)

#.................................................
#.................................................
# Model scenarios ####

null <- crossvalidation(G ~ empty_model,
                        data = mydata,
                        k = k,
                        folds = folds,
                        alpha = alpha)

null

env <- crossvalidation(G ~ minNT_veg + minNT_sow2rep + MLDS_rep,
                       data = mydata,
                       k = k,
                       folds = folds,
                       minsize = minsize,
                       alpha = alpha)


gen <- crossvalidation(G ~ minNT_veg + minNT_sow2rep,
                       data = mydata,
                       k = k,
                       folds = folds,
                       minsize = minsize,
                       alpha = alpha,
                       normal = prior,
                       gamma = gamma)


loc <- crossvalidation(G ~ lon + lat + xy + yx,
                       data = mydata,
                       k = k,
                       folds = folds,
                       minsize = minsize,
                       alpha = alpha)

loc

save(null, env, gen, loc,
     file = paste0(output, "models.rda"))


