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

#set.seed(123)
#k <- n
#folds <- sample(rep(1:k, times = ceiling(n / k), length.out = n))
folds <- as.integer(as.factor(ind$year))
k <- max(folds)

keep <- grepl("DT_|NT_|MLWS_|MLDS_", names(ind))

data <- cbind(G, empty_model = 0, ind[keep])

data <- cbind(G, ind[keep])

data <- data[,c(1:6)]
formula <- G ~ .
dots <- list(npseudo = 10,
             gamma = TRUE)
select.by = NULL
akaike.weights = TRUE
ncores = NULL
packages = NULL

# f <- forward(G ~ .,
#              data = data, 
#              k = k,
#              folds = folds,
#              npseudo = 5,
#              gamma = TRUE)
# 
# 
# 
# # cvnull <- crossvalidation(G ~ empty_model,
# #                           data = data,
# #                           k = k,
# #                           folds = folds,
# #                           npseudo = npseudo,
# #                           alpha = alpha
# #                           
# # args <- list(data = data, 
# #              k = k, 
# #              folds = folds,
# #              minsize = minsize,
# #              gamma = gamma,
# #              alpha = alpha,
# #              npseudo = npseudo)
# # 
# # # create cluster to do parallelisation
# # ncores <- detectCores() -1
# # cluster <- parallel::makeCluster(ncores)
# # doParallel::registerDoParallel(cluster)
# # 
# # vars <- names(data)[3:ncol(data)]
# # fs <- length(vars)
# # i <- 1:fs
# # 
# # # get predictions from nodes and put in matrix
# # models <- try(foreach::foreach(i = i,
# #                                .combine = .comb) %dopar% (.forward_dopar(as.formula(
# #                                  paste0("G ~ ", paste(vars[i], collapse = " + "))
# #                                ),
# #                                args)))
# # 
# # 
# # write.csv(models, paste0(output, "model_estimates.csv"))
# # 
# # # Stop cluster connection
# # parallel::stopCluster(cluster)
# 
# 
# #.................................................
# #.................................................
# # Model scenarios ####
# cvnull <- crossvalidation(G ~ empty_model,
#                           data = data,
#                           k = k,
#                           folds = folds,
#                           npseudo = npseudo,
#                           alpha = alpha)
# 
# cvclim <- crossvalidation(G ~ minNT_sow2gra + minDT_veg,
#                           data = data,
#                           k = k,
#                           folds = folds,
#                           gamma = gamma,
#                           minsize = minsize,
#                           npseudo = npseudo,
#                           alpha = alpha)
# 
# cvclimgen <- crossvalidation(G ~ minNT_sow2gra + minDT_veg,
#                           data = data,
#                           k = k,
#                           folds = folds,
#                           gamma = gamma,
#                           normal = prior,
#                           minsize = minsize,
#                           npseudo = npseudo,
#                           alpha = alpha)
# 
# 
# keep <- grepl("DT_|NT_|DTR|year|lon|lat", names(ind))
# 
# data <- cbind(G, empty_model = 0, ind[keep])
# 
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
# capture.output(cat("Null model \n"),
#                formula(null),
#                pseudoR2(null),
#                cat("\n\n\n"),
#                cat("Design \n"),
#                formula(design),
#                pseudoR2(design),
#                cat("\n\n\n"),
#                cat("Climate \n"),
#                formula(clim),
#                pseudoR2(clim),
#                cat("\n\n\n"),
#                cat("Climate + Genes \n"),
#                formula(climgen),
#                pseudoR2(climgen),
#                file = paste0(output, "model_estimates.txt"))
# 
# # predict(null)
# # predict(design)
# # predict(clim)
# # predict(climgen)