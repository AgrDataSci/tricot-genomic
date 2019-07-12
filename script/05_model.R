# 
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
# remove NAs in grain yield
df %>% 
  filter(!is.na(gy_gm)) -> 
  gy

# only keep strict rankings of at least 2 distinct items
gy %>%
  group_by(id) %>%
  summarise(keep = length(id)) %>%
  mutate(keep = keep > 1) %>% 
  filter(keep) ->
  keep

# apply the logical vector
id <- gy$id %in% keep$id

# keep selected observations
df <- gy[id,]

id <- ind$id %in% keep$id

ind <- ind[id, ]

G <- to_rankings(data = df,
                 items = "genotype",
                 input = "gy_gm",
                 id = "id", 
                 grouped.rankings = TRUE)

P <- to_paircomp(G)

w <- G[1:length(G) , , as.grouped_rankings = FALSE]
w[w==0] <- NA
w <- apply(w, 1, function(x) sum(!is.na(x)))
w <- 1 + (w * 0.1)


mod <- PlackettLuce(G)

# normal prior object from additive matrix
prior <- list(mu = coef(mod),
              Sigma = additivemat)


# mod2 <- PlackettLuce(G, normal = prior, gamma = TRUE)
# 
# AIC(mod)
# AIC(mod2)


#.................................................
#.................................................
# Set parameters for forward selection ####
out <- !grepl("Rx1day|Rx5day|SU|Rtotal|SDII", names(ind))
ind <- ind[out]

data <- cbind(G, empty_model = 0, ind[-ncol(ind)])

n <- nrow(data)

# names(data)
# 
# head(data)
# 
# 
# null <- pltree(G ~ empty_model,
#                data = data, 
#                weights = w)
# 
# AIC(null)
# 
# 
# tree <- pltree(G ~ minNT_rep + minNT_sow2rep + minNT_sow2gra, 
#                data = data, 
#                weights = w,
#                normal = prior, 
#                gamma = TRUE,
#                bonferroni = TRUE)
# 
# tree
# 
# AIC(null)
# AIC(tree)
# 
# pseudoR2(null, newdata = data)
# 
# pseudoR2(tree, newdata = data)
# 
# logLik(tree)

# 
# AIC(tree)
# 
# pseudoR2(tree, newdata = data)
# 
# null <- pltree(G ~ empty_model, data = data, weights = w, gamma = TRUE, normal = prior)
# 
# AIC(null)
# 
# pseudoR2(tree, newdata = data)
# pseudoR2(null, newdata = data)

# cross-validation parameters
set.seed(123)
k <- 100
folds <- sample(rep(1:k, times = ceiling(n / k), length.out = n))


library(psychotree)

null <- crossvalidation(G ~ empty_model, 
                        data = data, 
                        k = k, 
                        folds = folds,
                        bonferroni = TRUE)

clim <- crossvalidation(G ~ minNT_rep + minNT_sow2rep + minNT_sow2gra,
                        data = data, 
                        k = k, 
                        folds = folds,
                        bonferroni = TRUE,
                        normal = prior, 
                        gamma = TRUE)

null
clim


# PlackettLuce parameters
minsize <- round((nrow(data)*0.15), -1)
bonferroni <- TRUE
alpha <- 0.05
normal <- prior
gamma <- TRUE
npseudo <- 0.5

# forward selection parameters 
vars <- names(data)[grepl("minNT|maxNT", names(data))] #names(data)[2:27]#ncol(data)]
vars <- c("empty_model", vars)
var_keep <- character()
coeffs <- list()
counter <- 1
best <- TRUE

# parallelisation parameters
ncores <- 2
cluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cluster)

# keep running until the model get its best performance
while (best) {
  
  cat("\nForward Selection. Step ", counter, "\n Time: ", date(), "\n")
  
  fs <- length(vars)
  
  args <- list(data = data, 
               k = k, 
               folds = folds, 
               alpha = alpha,
               minsize = 250, 
               npseudo = 20,
               bonferroni = TRUE)
  
  i <- 1:fs
  
  # get predictions from nodes and put in matrix
  models <- foreach::foreach(i = i,
                             .combine = .comb) %dopar% (.forward_dopar(as.formula(
                               paste0("G ~ ", paste(c(var_keep, vars[i]), collapse = " + "))
                                 ),
                                 args))
  
  dimnames(models) <- list(1:fs,
                           paste0("bin",1:k), 
                           c("AIC","deviance","logLik","MaxLik","CraggUhler", "Agresti"))
  
  # take the matrix with selected goodness of fit
  modpar <- models[, , dimnames(models)[[3]] %in% "deviance"]
  
  if(is.null(dim(modpar))) {
    modpar <- t(as.matrix(modpar))
  }
  
  
  # if TRUE
  # then the highest value is taken 
  if (aw) {
    
    # calculate akaike weigths 
    # adjust to function to the matrix arrangement
    if (nrow(modpar) > 1) {
      modpar <- apply(modpar, 2, function(x) {
        akaike_weights(x)[[3]]
      })
    } else {
      modpar <- apply(modpar, 1, function(x) {
        akaike_weights(x)[[3]]
      })
      # turn it into a matrix again
      modpar <- t(as.matrix(modpar))
    }
    
    # then take the foldsize mean
    modpar <- apply(modpar, 1, function(x) {
      .mean_crossvalidation(object = x, 
                            folds = folds,
                            mean.method = "foldsize")
    })
    
    index_best <- which.max(modpar)
    
    value_best <- modpar[index_best]
    
    best <- .is_greater(value_best, baseline)
    
    
  }
  
  # if FALSE
  # select accordingly to the chosen method
  if (!aw) {
    
    modpar <- apply(modpar, 1, function(x){
      .mean_crossvalidation(x, 
                            folds = folds, ...)
    })
    
    # if AIC or deviance are selected then the model 
    # with lower value is the best 
    # other methods take the higher value
    if (select.by %in% c("AIC","deviance")) {
      
      index_best <- which.min(modpar)
      
      value_best <- modpar[index_best]
      
      best <- .is_lower(value_best, baseline)
      
    } else {
      
      index_best <- which.max(modpar)
      
      value_best <- modpar[index_best]
      
      best <- .is_greater(value_best, baseline)
    }
  }
  
  # refresh baseline 
  baseline <- value_best
  
  # take the name of best variable
  best_model <- exp_var[index_best]
  
  if (best_model == "empty_model") { 
    best <- FALSE 
    var_keep <- "empty_model"
  }
  
  # this is to prevent error in the array dimension 
  # if all variables are included in the model
  if (length(dim(models)) == 3) {
    models_avg <- apply(models, c(1,3), function(x){
      .mean_crossvalidation(x, 
                            folds = folds)
    })
    
  }else{
    models_avg <- apply(models, 2, function(x){
      .mean_crossvalidation(x, 
                            folds = folds)
    })
  }
  
  # model calls to add into list of parameters
  call_m <- paste0(Y, " ~ ", paste(paste(var_keep, collapse = " "), exp_var))
  call_m <- tibble::as_tibble(cbind(call = call_m, models_avg))
  call_m[2:ncol(call_m)] <- lapply(call_m[2:ncol(call_m)], as.numeric)
  call_m <- list(call_m, models)
  
  # take outputs from this run and add it to the list of parameters
  coeffs[[counter]] <- call_m
  
  if (best) {
    
    # remove best variable for the next run
    exp_var <- exp_var[!grepl(best_model, exp_var)]
    
    # remove empty var from the first run, no longer necessary
    exp_var <- exp_var[!grepl("empty_model", exp_var)]
    
    # keep this model for the next run
    var_keep <- c(var_keep, best_model)
    
    cat(" Best model found:", 
        paste0(Y, " ~ ", paste(var_keep, collapse = " + ")), 
        "\n\n")
    
  }
  
  # update counter (number of runs in 'while')
  counter <- counter + 1
  
  # prevent while loop to broke if the model fits with all variables
  if(length(exp_var) == 0) {
    best <- FALSE
  }
  
}

# Stop cluster connection
parallel::stopCluster(cluster)
