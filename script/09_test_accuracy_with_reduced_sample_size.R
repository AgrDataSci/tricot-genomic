# Model scenarios with explanatory variables 
# to predict variety performance using PlackettLuce models
library("tidyverse")
library("magrittr")
library("PlackettLuce")
library("gosset")
library("patchwork")
library("caret")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
                file = "script/session_info/09_test_accuracy_with_reduced_sample_size_session_info.txt")

# this function will sample the data in a balanced way
# for example when a group has much more observation than other 
sample2 <- function(x, f, ratio = NULL, seed = NULL) {
  
  x <- cbind(x, f)
  
  x <- split(x[,1], f = x[,2])
  
  set.seed(seed)
  x <- lapply(x, function(y){
    sample(y, 
           size = (length(y) * ratio), 
           replace = FALSE)
  })
  
  x <- as.vector(unlist(x))
  
  return(x)
  
}



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
                  group = TRUE)


R <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id",
                  group = FALSE)


#.................................................
#.................................................
# normal prior coefficients ####
mod <- PlackettLuce(G)

items <- names(coef(mod))

genitems <- items[items %in% dimnames(gene)[[1]]]

keep <- dimnames(gene)[[1]] %in% genitems

# keep only the items used in this analysis
gene <- gene[keep, ]

# remove those with no variance
drop <- nearZeroVar(gene)
gene <- gene[, -drop]

set.seed(765)
am <- matrix(runif(length(items) * ncol(gene), -1e-1, 1e-1), 
             nrow = length(items), 
             ncol = dim(gene)[[2]],
             dimnames = list(items,
                             dimnames(gene)[[2]]))

am[genitems,] <- gene

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
bonferroni <- FALSE
gamma <- TRUE
folds <- as.numeric(as.factor(ind$year))
k <- max(folds)


# select variables
keep <- which(grepl("DT_|NT_|MLDS_rep|lon|lat|xy|yx|year", names(ind)))

ind <- ind[, keep]

# put all together as mydata
mydata <- cbind(G, ind)


#.................................................
#.................................................
# Model with OA ####
gen <- crossvalidation(G ~ minNT_veg + maxNT_rep,
                       data = mydata,
                       k = k,
                       folds = folds,
                       minsize = minsize,
                       alpha = alpha,
                       normal = prior,
                       gamma = gamma)

gen$raw$estimators


# make scenarios with a small part of the data
# what happens with we take 10, 25, 50 and 75% of the data
# first we create a vector that is the seeds for the 100 runs in 
# each scenario
set.seed(10387)
s <- as.integer(runif(500, 2, 10000))
out <- rep(c(0.25, 0.5, 0.75, 0.85, 0.95), times = 100)

# then we run over these seeds to sample the data and take the 
# Kendall tau 
kt <- data.frame()

for(i in seq_along(s)) {
  print(i)
  si <- sample2(rownames(mydata), mydata$year, 1 - out[i], seed = s[i])
  dt <- mydata[si,]
  minsize <- round(nrow(dt) * 0.2, -2)
  k <- length(unique(dt$year))
  folds <- as.numeric(as.factor(dt$year))
  
  cv <- crossvalidation(G ~ minNT_veg + maxNT_rep,
                        data = dt,
                        k = k,
                        folds = folds,
                        minsize = minsize,
                        alpha = alpha,
                        normal = prior,
                        gamma = gamma)
  
  cv$coeffs$sample <- out[i]
  
  kt <- rbind(kt, cv$coeffs[, c("sample", "kendallTau")])
  
}

kt <- kt[order(kt$sample), ]

kt$run <- rep(1:100, nrow(kt)/100)

kt$trait <- "OA"

write_csv(kt, "output/reduce_sample_size/kendall_tau_smaller_samples.csv")

kt %>% 
  group_by(sample) %>% 
  summarise(kendall = mean(kendallTau))



# Now using GY
#...........................
# grain yield from decentralized
# remove NAs in grain yield
gy <- rank_numeric(df,
                   "genotype",
                   "gy_gm",
                   "id")

keep <- !is.na(gy)

gy <- gy[keep, ]

gy_ind <- ind[keep, ]

GY <- group(gy, index = 1:length(gy))

dataGY <- cbind(GY, gy_ind)

folds_gy <- folds[keep]

gen_gy <- crossvalidation(GY ~ minNT_veg + maxNT_rep,
                          data = dataGY,
                          k = k,
                          folds = folds_gy,
                          minsize = minsize,
                          alpha = alpha,
                          normal = prior)
gen_gy

# then we run over these seeds to sample the data and take the 
# Kendall tau 
ktgy <- data.frame()

for(i in seq_along(s)) {
  print(i)
  si <- sample2(rownames(dataGY), dataGY$year, 1 - out[i], seed = s[i])
  dt <- dataGY[si, ]
  minsize <- round(nrow(dt) * 0.3, -2)
  k <- length(unique(dt$year))
  folds <- as.numeric(as.factor(dt$year))
  
  cv <- crossvalidation(GY ~ minNT_veg + maxNT_rep,
                        data = dt,
                        k = k,
                        folds = folds,
                        minsize = minsize,
                        alpha = alpha,
                        normal = prior)
  
  cv$coeffs$sample <- out[i]
  
  ktgy <- rbind(ktgy, cv$coeffs[, c("sample", "kendallTau")])
  
}

ktgy <- ktgy[order(ktgy$sample), ]

ktgy$run <- rep(1:100, nrow(ktgy)/100)

ktgy$trait <- "GY"

ktgy

write_csv(ktgy, "output/reduce_sample_size/kendall_tau_smaller_samples_GY.csv")

kendall <- rbind(kt, ktgy)

kendall$n_plots <- floor((nrow(mydata)*4/41) * (1- kendall$sample))

write_csv(kendall, "output/reduce_sample_size/kendall_tau_OA_GY.csv")

kendall %>% 
  group_by(trait, sample, n_plots) %>% 
  summarise(kendall = mean(kendallTau), 
            sd = sd(kendallTau)) %>% 
  ungroup()






