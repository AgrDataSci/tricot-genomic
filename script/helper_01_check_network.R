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


# not genotyped varieties
out <- "data/raw/not_in_genotyping.csv"
out %<>% 
  read_csv() %>% 
  select(accession)

vars <- bind_cols(accession = unlist(df[paste0("variety_",letters[1:4])]),
                  region = rep(df$region, 4))


vars %<>%
  group_by(region) %>% 
  distinct(accession)


vars %<>% 
  bind_cols(. , out = vars$accession %in% out$accession)

vars %>%  
  filter(out == TRUE) ->
  out


vars %>% 
  filter(out == FALSE) ->
  inn


R <- to_rankings(items = df[paste0("variety_",letters[1:4])],
                 input = df[paste0("rank_variety_",letters[1:4])])


R <- R[, dimnames(R)[[2]] %in% inn$accession]


adj <- adjacency(R)

adj <- as.vector(adj)

adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))

dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])

adj <- btdata(adj, return_graph = TRUE)

plot.igraph(adj$graph, vertex.size = 10, edge.arrow.size = 0.1) 


mod <- PlackettLuce(R)
summary(mod)
plot(qvcalc(mod), las = 2, cex.axis = 0.7)
