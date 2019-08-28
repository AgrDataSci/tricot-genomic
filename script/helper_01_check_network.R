library("tidyverse")
library("magrittr")
library("gosset")
library("PlackettLuce")
library("BradleyTerryScalable")
library("igraph")

# full dataset
df <- "data/durumwheat.csv"
df %<>%
  read_csv()


table(df$accession, df$region)

R <- rank_PL(data = df,
             items = "genotype",
             input = "farmer_rank",
             id = "id")


adj <- adjacency(R)

connectivity(adj)

adj <- as.vector(adj)

adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))

dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])

adj <- btdata(adj, return_graph = TRUE)


svg(filename = "output/network.svg", width = 12, height = 12, pointsize = 20)
plot.igraph(adj$graph, vertex.size = 5, edge.arrow.size = 0.1) 
dev.off()

mod <- PlackettLuce(R)
summary(mod)
plot(qvcalc(mod), las = 2, cex.axis = 0.7)
