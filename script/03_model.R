library(gosset)
library(PlackettLuce)
library(BradleyTerryScalable)
library(igraph)

R <- to_rankings(items = df[paste0("variety_",letters[1:4])],
                 input = df[paste0("rank_variety_",letters[1:4])])


adj <- PlackettLuce::adjacency(R)

adj <- as.vector(adj)

adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))

dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])

adj <- btdata(adj, return_graph = TRUE)

plot.igraph(adj$graph, vertex.size = 10, edge.arrow.size = 0.1) 


mod <- PlackettLuce(R)
summary(mod)
plot(qvcalc(mod), las = 2, cex.axis = 0.7)
