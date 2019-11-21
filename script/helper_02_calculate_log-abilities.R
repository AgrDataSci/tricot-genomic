# Get log abilities from PlackettLuce and Bradley-Terry models 
# as a proxy for genomic prediction

library("tidyverse")
library("magrittr")
library("gosset")
library("PlackettLuce")
library("BradleyTerry2")
library("svglite")

#.....................................
#.....................................
# Read data ####

df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

output <- "output/log-abilities/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

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
# probabilities for all years ####

R <- rank_PL(data = df,
             items = "genotype",
             input = "farmer_rank",
             id = "id")

mod <- PlackettLuce(R)

pars <- itempar(mod, log = FALSE)

pars <- bind_cols(genotype = names(pars),
                  pow_rank = as.vector(pars),
                  season = rep("all", length(pars)))


#.....................................
#.....................................
# probabilities for each year ####

# split into seasons
R <- split(df, df$year)

probs <- array(NA, dim = c(nrow(pars), 3, 3),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_rank","season"),
                               1:3))
# run over years
for(i in seq_along(R)) {
  r <- rank_PL(data = R[[i]],
               items = "genotype",
               input = "farmer_rank",
               id = "id")
  
  mod <- PlackettLuce(r)
  
  p <- itempar(mod, log = FALSE)
  
  p <- as.matrix(p)
  
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  
  probs[,1,i] <- pars$genotype
  
  probs[, 3, i] <- names(R[i]) 
  

}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
    mutate(pow_rank = as.numeric(pow_rank))

probs <- bind_rows(pars, probs)

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_rank = ifelse(is.na(pow_rank), mean(pow_rank, na.rm = TRUE), pow_rank)) %>% 
  ungroup()

probs %<>% 
  group_by(season) %>% 
  mutate(pow_rank = pow_rank / sum(pow_rank)) %>% 
  ungroup()


winprobs <- probs


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

# ............................
# ............................
# all data 
R <- rank_PL(data = gy,
             items = "genotype",
             input = "gy_gm",
             id = "id")

mod <- PlackettLuce(R)

pars <- itempar(mod, log = FALSE)

pars <- bind_cols(genotype = names(pars),
                  pow_rank = as.vector(pars),
                  season = rep("all", length(pars)))



# .....................................
#.....................................
# probabilities for each year ####

# split into seasons
R <- split(gy, gy$year)

probs <- array(NA, dim = c(nrow(pars), 3, 3),
               dimnames = list(pars$genotype, 
                               c("genotype","pow_rank","season"),
                               1:3))
# run over years
for(i in seq_along(R)) {
  r <- rank_PL(data = R[[i]],
               items = "genotype",
               input = "gy_gm",
               id = "id")
  
  mod <- PlackettLuce(r)
  
  p <- itempar(mod, log = FALSE)
  
  p <- as.matrix(p)
  
  probs[dimnames(p)[[1]], 2, i] <- p[,1]
  
  probs[,1,i] <- pars$genotype
  
  probs[, 3, i] <- names(R[i]) 
  
  
}

probs <- as_tibble(apply(probs, 2L, c))

probs %<>% 
  mutate(pow_rank = as.numeric(pow_rank))

probs <- bind_rows(pars, probs)

probs <- probs[probs$genotype %in% items, ]

probs %<>% 
  group_by(genotype) %>% 
  mutate(pow_rank = ifelse(is.na(pow_rank), mean(pow_rank, na.rm = TRUE), pow_rank)) %>% 
  ungroup()

probs %<>% 
  group_by(season) %>% 
  mutate(pow_rank = pow_rank / sum(pow_rank)) %>% 
  ungroup()

probs %<>% 
  rename(pow_gy = pow_rank)

winprobs %<>% 
  select(pow_rank) %>% 
  bind_cols(. , probs) %>% 
  select(genotype, season, pow_rank, pow_gy)



write_csv(winprobs, paste0(output, "log-abilities.csv"))



# # ...................................
# # ...................................
# # farmer rank vs yield ####
# 
# # yield rankings into a parsed matrix
# YR <- YR[1:length(YR),,as.grouped_rankings = FALSE]
# 
# # farmer ranking into a PL object
# FR <- to_rankings(data = gy,
#                   items = "genotype",
#                   input = "farmer_rank",
#                   id = "id", 
#                   grouped.rankings = TRUE)
# 
# # then into a parsed matrix 
# FR <- FR[1:length(FR), , as.grouped_rankings = FALSE]
# 
# kendall <- kendallTau(FR, YR)
# 
# write_csv(kendall, paste0(output, "kendall_correlation.csv"))
# dimnames(R)[2]
# 
# adj <- PlackettLuce::adjacency(R)
# 
# adj <- as.vector(adj)
# 
# adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))
# 
# dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])
# 
# adj <- btdata(adj, return_graph = TRUE)
# 
# plot.igraph(adj$graph, vertex.size = 10, edge.arrow.size = 0.1)
# 
# 
# mod <- PlackettLuce(R)
# summary(mod)
# plot(qvcalc(mod), las = 2, cex.axis = 0.7)

