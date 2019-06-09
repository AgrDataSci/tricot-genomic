library("tidyverse")
library("magrittr")
library("gosset")
library("PlackettLuce")


df <- "data/durumwheat.csv"
df %<>% 
  read_csv()


output <- "output/exploring/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

#.....................................
#.....................................
# farmer rank ####

# create grouped rankings 
G <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)

# add explanatory vars
# to_rankings internaly order the rankings by their id 
# so I will keep the unique values by id and combine with the 
# grouped_rankings

df %>% 
  arrange(id) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  select(year, region) -> 
  expvar

G <- cbind(G, expvar)

mod <- pltree(G ~ ., data = G)

winprobs <- itempar(mod, log = TRUE)

winprobs <- bind_cols(genotype = dimnames(winprobs)[[2]],
                      pow_rank = as.vector(winprobs))


#...........................
#...........................
# grain yield ####

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
gy <- gy[id,]


YR <- to_rankings(data = gy,
                  items = "genotype",
                  input = "gy_gm",
                  id = "id", 
                  grouped.rankings = TRUE)


mod_gy <- pltree(YR ~ 1, data = YR)


winprobs_gy <- itempar(mod_gy, log = TRUE)

winprobs_gy <- bind_cols(genotype = dimnames(winprobs_gy)[[2]],
                         pow_gy = as.vector(winprobs_gy))

winprobs %<>% 
  merge(. , winprobs_gy, all.x = TRUE) %>% 
  as_tibble()


write_csv(winprobs, paste0(output, "probability_of_winning.csv"))


# ...................................
# ...................................
# farmer rank vs yield ####

# yield rankings into a parsed matrix
YR <- YR[1:length(YR),,as.grouped_rankings = FALSE]

# farmer ranking into a PL object
FR <- to_rankings(data = gy,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id", 
                  grouped.rankings = TRUE)

# then into a parsed matrix 
FR <- FR[1:length(FR), , as.grouped_rankings = FALSE]

kendall <- kendallTau(FR, YR)

write_csv(kendall, paste0(output, "kendall_correlation.csv"))


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

