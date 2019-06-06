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
                 items = "accession",
                 input = "farmer_rank",
                 id = "id", 
                 grouped.rankings = TRUE)



mod <- pltree(R ~ 1, data = R, bonferroni = TRUE)

summary(mod)

plot_nodes(mod)




df %>% 


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
