######################
# this script calculates logabilities from farmer fields
# Date: 01/05/2020
######################

library(RCurl)
library(tidyr)
library(dplyr)
library(gosset)
library(PlackettLuce)
library(BradleyTerry2)
library(svglite)
library(magrittr)
library(readr)


# Read data ####
#load the original dataset
dfurl<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/durumwheat.csv")
df<-read.csv(text=dfurl)

#create output folder
output <- "output/log-abilities/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

#.....................................
#.....................................
# Filter data ####
# keep only genotyped assessions
items <- unique(df$genotype)

items <- items[grepl("_D", items)]

df <- df[df$genotype %in% items, ]

# only keep strict rankings of at least 2 distinct items
df %>%
  group_by(id) %>%
  summarise(keep = length(id)) %>%
  mutate(keep = keep > 1) %>% 
  filter(keep) ->
  keep

# apply the logical vector
id <- df$id %in% keep$id

# keep selected observations
df <- df[id,]

#.....................................
#.....................................
# Add grouping variables ####

#load environmental quantiles derived from other script
load("output/environmental.quantiles.by.centralized.loc.Rdata")
head(envqt)

#attach reference to locations
envqt[,2]<-sub("^", "ger.", envqt[,2])
envqt[,3]<-sub("^", "hs.", envqt[,3])

#include them in the overall dataset
df<-merge(df, envqt, by.x="id", by.y="id", all.x=T)

#.....................................
#.....................................
# Compute probabilities for all years ####
# on overall rank

R <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id",
                  ascending = TRUE)

mod <- PlackettLuce(R)

pars <- coef(mod, log = FALSE)

pars <- bind_cols(genotype = names(pars),
                  pow_rank = as.vector(pars),
                  group = rep("all", length(pars)))

#.....................................
#.....................................
# probabilities for each year ####

# split into seasons
R <- split(df, df$year)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_rank","group"),
                               1:length(R)))
# run over years
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "farmer_rank", 
                    id = "id",
                    ascending = TRUE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
  
  
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_rank = as.numeric(pow_rank))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_rank = ifelse(is.na(pow_rank), mean(pow_rank, na.rm = TRUE), pow_rank)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_rank = pow_rank / sum(pow_rank)) %>% 
  ungroup()

winprobs1 <- probs

#.....................................
#.....................................
# probabilities for geregera environmental quantile ####

# split into quantiles for geregera
R <- split(df, df$ger.qt)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_rank","group"),
                               1:length(R)))
# run over quantiles
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "farmer_rank", 
                    id = "id",
                    ascending = TRUE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_rank = as.numeric(pow_rank))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_rank = ifelse(is.na(pow_rank), mean(pow_rank, na.rm = TRUE), pow_rank)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_rank = pow_rank / sum(pow_rank)) %>% 
  ungroup()

winprobs2 <- probs


#.....................................
#.....................................
# probabilities for hagreselam quantiles ####

# split into quantiles for hagreselam
R <- split(df, df$hs.qt)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_rank","group"),
                               1:length(R)))
# run over quantiles
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "farmer_rank", 
                    id = "id",
                    ascending = TRUE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_rank = as.numeric(pow_rank))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_rank = ifelse(is.na(pow_rank), mean(pow_rank, na.rm = TRUE), pow_rank)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_rank = pow_rank / sum(pow_rank)) %>% 
  ungroup()

winprobs3 <- probs


#.....................................
#.....................................
# put everything together and write it on file####

winprobs <- bind_rows(pars, winprobs1,winprobs2, winprobs3)

out <- winprobs %>% spread(group, pow_rank)

write_csv(out, paste0(output, "pow.log-abilities.csv"))

#############################################
#     DO THE SAME WITH GY
#############################################
#...........................
#...........................
# grain yield ####

# remove NAs in grain yield
df %>% 
  filter(!is.na(gy_gm)) -> 
  gy

gy %<>% 
  group_by(id) %>% 
  distinct(gy_gm, .keep_all = TRUE) %>% 
  ungroup()

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
gy <- gy[id,]


#.....................................
#.....................................
# Compute probabilities for all years ####
# on gy rank

R <- rank_numeric(data = gy,
                  items = "genotype",
                  input = "gy_gm",
                  id = "id",
                  ascending = FALSE)

mod <- PlackettLuce(R)

gy.pars <- coef(mod, log = FALSE)

gy.pars <- bind_cols(genotype = names(gy.pars),
                  pow_gy = as.vector(gy.pars),
                  group = rep("all", length(gy.pars)))

#.....................................
#.....................................
# probabilities for each year ####

# split into seasons
R <- split(gy, gy$year)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_gy","group"),
                               1:length(R)))
# run over years
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "gy_gm",
                    id = "id", 
                    ascending = FALSE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
  
  
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_gy = as.numeric(pow_gy))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_gy = ifelse(is.na(pow_gy), mean(pow_gy, na.rm = TRUE), pow_gy)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_gy = pow_gy / sum(pow_gy)) %>% 
  ungroup()

gy.winprobs1 <- probs

#.....................................
#.....................................
# probabilities for geregera environmental quantile ####

# split into quantiles for geregera
R <- split(gy, gy$ger.qt)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_gy","group"),
                               1:length(R)))
# run over quantiles
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "gy_gm",
                    id = "id", 
                    ascending = FALSE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_gy = as.numeric(pow_gy))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_gy = ifelse(is.na(pow_gy), mean(pow_gy, na.rm = TRUE), pow_gy)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_gy = pow_gy / sum(pow_gy)) %>% 
  ungroup()

gy.winprobs2 <- probs


#.....................................
#.....................................
# probabilities for hagreselam quantiles ####

# split into quantiles for hagreselam
R <- split(gy, gy$hs.qt)

probs <- array(NA, dim = c(nrow(pars), 3, length(R)),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_gy","group"),
                               1:length(R)))
# run over quantiles
for(i in seq_along(R)) {
  r <- rank_numeric(data = R[[i]],
                    items = "genotype",
                    input = "gy_gm",
                    id = "id", 
                    ascending = FALSE)
  mod <- PlackettLuce(r)
  p <- coef(mod, log = FALSE)
  p <- as.matrix(p)
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  probs[,1,i] <- pars$genotype
  probs[, 3, i] <- names(R[i]) 
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_gy = as.numeric(pow_gy))

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_gy = ifelse(is.na(pow_gy), mean(pow_gy, na.rm = TRUE), pow_gy)) %>% 
  ungroup()

probs %<>% 
  group_by(group) %>% 
  mutate(pow_gy = pow_gy / sum(pow_gy)) %>% 
  ungroup()

gy.winprobs3 <- probs

# put everything together and write it on file####

gy.winprobs <- bind_rows(gy.pars, gy.winprobs1,gy.winprobs2, gy.winprobs3)

gy.out <- gy.winprobs %>% spread(group, pow_gy)

write_csv(out, paste0(output, "gy.log-abilities.csv"))

############################
#merge both datasets into one
colnames(out)<-sub("$", ".pow", colnames(out))
colnames(gy.out)<-sub("$", ".gy", colnames(gy.out))

logab<-merge(out, gy.out, by.x="genotype.pow", by.y="genotype.gy")
rownames(logab)<-logab[,1]
logab<-logab[,-1]

save(logab, file=paste0(output, "all.traits.log-abilities.Rdata"))
