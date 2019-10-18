# compute gain in grain yield based on model predictions
reliability_gy <- function(object, gy, 
                           weights = c(0.5, 0.3, 0.2), 
                           sort = "yield") {
  
  if(sum(weights) != 1) {
    stop("Sum of weights is not 1")
  }
  
  items <- dimnames(object)[[2]]
  
  i <- gy$genotype
  
  object <- object[, items %in%  i ]
  
  # put it into rankings
  object <- t(apply(object, 1, function(X){
    t(gosset::rank_decimal(X)["rank"])
  }))
  
  dimnames(object)[[2]] <- i
  
  # identify the best three
  yield <- apply(object, 1, function(x) {
    # take the best three
    index <- which(x %in% c(1:3))
    # take the names of the best three
    best <- x[index]
    
    # get the yield for the best three
    if (sort=="yield") {
      bestnames <- names(best)
      y <- gy[gy$genotype %in% bestnames, "gy" ]
      y <- sort(t(y), decreasing = TRUE)
    }
    if (sort=="rank") {
      bestnames <- names(sort(best))
      y <- gy[gy$genotype %in% bestnames, "gy" ]
    }
    
    # weight it
    y <- y * weights
    # sum it
    sum(y)
    
  })
  
  return(yield)
  
}


# Compute reliability, relative gains compared to a reference object
# return a raster object with relative gains ranging from 0 to 1
reliability <- function(object, reference = NULL, 
                        alternative = NULL, 
                        xy = NULL, 
                        avg.probs = NULL, 
                        resolution = c(0.04166667, 0.04166667), 
                        projection = "+proj=longlat +datum=WGS84",
                        convexhull = TRUE){
  
  #if average probabilities from cross-validation is not provided then use a simple matrix from "predict" function
  if(is.null(avg.probs))  {winProbs <- predict(object)  }
  if(!is.null(avg.probs)) {winProbs <- avg.probs}
  
  #get number os rows and item names
  n <- dim(winProbs)[1]
  items <- dimnames(winProbs)[[2]]
  
  #check if the provided reference item exists within the item names
  if(!reference %in% items) {stop(paste0("unkown reference item, choose one among: ", toString(items), "\n"))}
  
  #if no alternative item is given, pick one using "worstRegret" function.
  if(is.null(alternative)) {
    WR <- worstRegret(tree)
    alternative <- WR[which.min(WR[,3]), 1]
  }
  #check if provided alternative item exists within the item names
  if(!alternative %in% items) stop(paste0("unkown alternative item, choose one among: ", toString(items), "\n"))
  
  #calculate the probability of the alternative item beat the reference item
  reliab <- winProbs[,alternative] /  (winProbs[,reference] + winProbs[,alternative])
  
  #Make the raster with reliability
  #geostats analysis
  #change names in given coordinates
  names(xy) <- c("x","y")
  myproj <- projection
  
  #Make a convex hull using coordinates
  #If dataset is larger than 5000 observations
  #take a sample to avoid problems in RAM memory
  if(n > 5000){
    #calculate largest distance between coordinates
    largedist <- xy[ sample(row.names(xy), 4000) , ]
    largedist <- raster::pointDistance(largedist, longlat = FALSE)
    largedist <- max(largedist, na.rm = TRUE)
    
  }
  if(n < 5000){
    #calculate largest distance between coordinates
    largedist <- raster::pointDistance(xy, longlat = FALSE)
    largedist <- max(largedist, na.rm = TRUE)
  }
  
  #convex hull
  hull <- dismo::convHull(xy, lonlat=TRUE)
  hull <- rgeos::gBuffer(hull@polygons, width=0.1*largedist)
  
  #set projection
  crs(hull) <- myproj
  #define extention for interpolation
  myext <- raster::extent(hull)
  
  #Add realiability values to cooordinates
  xy$reliab <- reliab
  
  #Define raster
  r <- raster::raster(hull)
  #Set the resolution of the cells
  raster::res(r) <- resolution
  
  #Resample coordinates and remove possible points within the same grid
  s <- dismo::gridSample(xy[,c("x","y")], r, n=1)
  xy <- xy[rownames(s),]
  
  #convert to coordinates
  sp::coordinates(xy) <- ~x+y
  sp::proj4string(xy) <- CRS(myproj)
  
  #Define the grid extent 
  x_range <- c(myext[1],myext[2])
  y_range <- c(myext[3],myext[4])
  
  #Create an empty grid where n is the total number of cells
  grd <- expand.grid(x = seq(x_range[1], x_range[2], by = 0.01),
                     y = seq(y_range[1], y_range[2], by = 0.01))
  
  sp::coordinates(grd) <- ~x+y
  sp::gridded(grd)     <- TRUE
  sp::proj4string(grd) <- proj4string(xy)
  
  #interpolate the grid cells
  p_idw <- gstat::idw(reliab ~ 1, locations = xy, newdata = grd)
  #convert to raster object then clip
  r <- raster::raster(p_idw)
  
  #Crop raster to fit the convexhull
  if(convexhull) {r <- raster::mask(r, hull)}
  
  r <- raster::crop(r, myext)
  
  return(r)
  
}

