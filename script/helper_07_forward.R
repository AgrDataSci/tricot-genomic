library("abind")
library("gosset")
library("parallel")
library("doParallel")

# combine results from parallel
.comb <- function(...) {
  abind::abind(..., along = 2, force.array = TRUE)
}

# model call for parallel
.forward_dopar <- function(formula, args){
  
  args <- c(formula, args)
  
  m <- do.call(gosset::crossvalidation, args)
  
  result <- m$raw$estimators$deviance

  return(result)
  
}


# logical function for > greater 
.is_greater <- function(x, y) {
  x > y
}


# logical function for < lower
.is_lower <- function(x, y) {
  x < y
}


# Compute weighted means in cross-validation
.mean_crossvalidation <- function(object, folds = NULL, 
                                  mean.method = NULL, 
                                  ...){
  # take length of folds
  N <- length(folds)
  
  if (is.null(mean.method)) {mean.method <- "stouffer"}
  
  # weight of imbalanced folds
  if (mean.method == "stouffer") {
    # take the number of folds
    max_folds <- max(folds)
    
    # make a table of folds and take
    # how many observations each fold has
    foldsize <- table(folds)
    
    # take the weight of each fold
    # first, the squared root of foldsize (observations per fold)
    # by the total number of observation
    wfold <- sqrt(as.vector(foldsize) / N)
    
    # then divide this vector by its sum
    wfold <- wfold / sum(wfold)
    # wfold <- (max_folds * wfold) / sum(wfold)
    
    # then we multiply the input values by the
    # weight of each fold
    stouffer <- object * wfold
    
    # sum these values and that is the stouffer mean
    mean <- sum(stouffer)
  }
  
  if (mean.method == "foldsize") {
    # make a table of folds and take
    # the number of observations per fold
    foldsize <- as.vector(table(folds))
    
    # fold size mean is the product of multiplication of object values by 
    # its number of observations then divided by the total number of observations
    mean <- sum(object * foldsize, na.rm = TRUE) / sum(foldsize, na.rm = TRUE)
  }
  
  if (mean.method == "equal") {
    mean <- mean(object, na.rm = TRUE)
  }
  
  return(mean)
  
}