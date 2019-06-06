library("tidyverse")
library("magrittr")
library("gosset")
library("PlackettLuce")
library("BradleyTerryScalable")
library("igraph")


df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

R <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)



mod <- pltree(R ~ 1, data = R, bonferroni = TRUE)

summary(mod)

plot_nodes(mod)

# check correlation between yield and farmer ranking
df %>% 
  filter(!is.na(gy_gm)) -> 
  yield

# only keep strict rankings of at least 2 distinct items
yield %>%
  group_by(farmer_no) %>%
  summarise(keep = n_distinct(accession) >= 2) %>%
  filter(keep) ->
  keep

# apply the logical vector
id <- yield$farmer_no %in% keep$farmer_no

# keep selected observations
yield <- yield[id,]


FR <- to_rankings(data = yield,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id", 
                  grouped.rankings = TRUE)

FR <- FR[1:length(FR), , as.grouped_rankings = FALSE]



YR <- to_rankings(data = yield,
                  items = "genotype",
                  input = "gy_gm",
                  id = "id", 
                  grouped.rankings = TRUE)

YR <- YR[1:length(YR), , as.grouped_rankings = FALSE]

kendallTau(FR, YR)




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

