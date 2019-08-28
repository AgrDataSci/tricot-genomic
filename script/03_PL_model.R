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
G <- rank_PL(data = df,
             items = "genotype",
             input = "farmer_rank",
             id = "id", 
             grouped.rankings = TRUE)

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

am[1:10, 1:10]

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
minsize <- 200 
npseudo <- 5
alpha <- 0.01
bonferroni <- TRUE
gamma <- TRUE
k <- 100
set.seed(123)
folds <- sample(rep(1:k, times = ceiling(n / k), length.out = n))
keep <- which(grepl("DT_|NT_|lon|lat|xy|yx", names(ind)))

data <- cbind(G, empty_model = rep(0, n), ind[keep])

#.................................................
#.................................................
# Model scenarios ####
out <- c(19, 24, 29, 49)

null <- crossvalidation(G ~ empty_model,
                        data = data,
                        drop.folds = out,
                        k = k,
                        folds = folds,
                        alpha = alpha)

print(null)

env <- crossvalidation(G ~ minDT_veg + minNT_veg,
                       data = data,
                       k = k,
                       folds = folds,
                       drop.folds = out,
                       minsize = minsize,
                       alpha = alpha)

print(env)

gen <- crossvalidation(G ~ minDT_veg + minNT_veg,
                       data = data,
                       k = k,
                       folds = folds,
                       drop.folds = out,
                       minsize = minsize,
                       alpha = alpha,
                       normal = prior)

print(gen)

loc <- crossvalidation(G ~ lon + lat + xy + yx,
                       data = data,
                       k = k,
                       folds = folds,
                       drop.folds = out,
                       minsize = minsize,
                       alpha = alpha)

print(loc)

save(null, env, gen, loc,
     file = paste0(output, "models.rda"))



