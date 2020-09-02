
library("analogues")
library("tidyverse")

dt <- read.csv("data/durumwheat.csv")

dt <- dt[!duplicated(dt$id), ]

dt <- dt[,c("region","lon","lat")]

st1 <- c(38.866667, 11.666667)

st2 <- c(39.166667 , 13.650000)

tmean <- getData("worldclim", var = "tmean", res = 2.5, lon = st2[1], lat = st2[2])
prec <- getData("worldclim", var = "prec", res = 2.5, lon = st2[1], lat = st2[2])
eth <- getData("GADM", country = "ETH", level = 1)

tmean <- stack(crop(tmean, eth))
prec <- stack(crop(prec, eth))

plot(tmean[[1]])
plot(eth, add = TRUE)
points(st1[1], st1[2], pch = "+")
points(st2[1], st2[2], pch = "+")
points(dt[,c("lon","lat")])


p <- list(st1, st2)
r <- list()

for(i in seq_along(p)){
  pars <- createParameters(x = p[[i]][1], 
                           y = p[[i]][2], 
                           vars = c("tmean","prec"),
                           weights = c(0.5,0.5),
                           ndivisions = c(12,12),
                           growing.season = c(5, 10),
                           rotation = "tmean",
                           threshold = 1,
                           env.data.ref = list(tmean, prec), 
                           env.data.targ = list(tmean, prec),
                           outfile="~/.",
                           fname=NA,
                           writefile=FALSE)
  
  sim <- calc_similarity(pars)
  
  r[[i]] <- sim
}

par(mfrow=c(2,1))
plot(r[[1]])
points(st1[1], st1[2], pch = "+")
points(dt[dt$region=="Amhara",c("lon","lat")])
plot(r[[2]])
points(st2[1], st2[2], pch = "+")
points(dt[dt$region=="Tigray",c("lon","lat")])


dt$s1 <- as.vector(raster::extract(r[[1]], dt[,c("lon","lat")]))
dt$s2 <- as.vector(raster::extract(r[[2]], dt[,c("lon","lat")]))

dt$s <- ifelse(dt$region == "Amhara", dt$s1,
               ifelse(dt$region == "Tigray", dt$s2,
                      ((dt$s1 + dt$s2) / 2)))


dt %>%
  group_by(region) %>%
  summarise(mean = mean(s),
            max = max(s),
            min = min(s))

boxplot(dt$s ~ dt$region)

# r <- stack(r)
# 
# r <- calc(r, mean)
# 
# dt$s <-  as.vector(raster::extract(r, dt[,c("lon","lat")]))
# 
# dt %>%
#   group_by(region) %>%
#   summarise(mean = mean(s),
#             max = max(s),
#             min = min(s))
# 
# boxplot(dt$s ~ dt$region)
