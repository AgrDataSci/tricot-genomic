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

#.....................................
#.....................................
# farmer rank ####

# create PlackettLuce rankings 
R <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id")

mod <- PlackettLuce(R)

svg(filename = paste0(output, "coeff_farmer_rank_PL.svg"),
    width = 10,
    height = 6.5,
    pointsize = 12)
plot(qvcalc(itempar(mod, log = TRUE)), las = 2)
dev.off()


winprobs <- itempar(mod, log = TRUE)

winprobs <- bind_cols(genotype = names(winprobs),
                      pow_rank = as.vector(winprobs))


#.....................................
#.....................................
# farmer rank BradleyTerry ####

# get a binomial rank
B <- rank_binomial(R, drop.null = TRUE)

mod <- BTm(cbind(win1, win2), player1, player2, ~ genotype,
           id = "genotype", data = B)


summary(mod)

x <- coef(mod)

names(x) <- gsub("genotype","",names(x))

x <- tibble(genotype = names(x),
            pow_rank_bt = x)

winprobs %<>% 
  merge(. , x, by = "genotype", all.x = TRUE) %>% 
  as_tibble()

winprobs %<>% 
  mutate(pow_rank_bt = ifelse(is.na(pow_rank_bt),
                              0,
                              pow_rank_bt))


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


# grain yield into rankings
YR <- to_rankings(data = gy,
                  items = "genotype",
                  input = "gy_gm",
                  id = "id")


mod_gy <- PlackettLuce(YR)


svg(filename = paste0(output, "coeff_grainyield_rank_PL.svg"),
    width = 10,
    height = 6.5,
    pointsize = 12)
plot(qvcalc(itempar(mod_gy, log = TRUE)), las = 2)
dev.off()


winprobs_gy <- itempar(mod_gy, log = TRUE)

winprobs_gy <- bind_cols(genotype = names(winprobs_gy),
                         pow_gy = as.vector(winprobs_gy))

winprobs %<>% 
  merge(. , winprobs_gy, all.x = TRUE) %>% 
  as_tibble()


#.....................................
#.....................................
# farmer rank BradleyTerry ####

# get a binomial rank
B <- rank_binomial(YR, drop.null = TRUE)

mod <- BTm(cbind(win1, win2), player1, player2, ~ genotype,
           id = "genotype", data = B)


summary(mod)

x <- coef(mod)

names(x) <- gsub("genotype","",names(x))

x <- tibble(genotype = names(x),
            pow_gy_bt = x)

winprobs %<>% 
  merge(. , x, by = "genotype", all.x = TRUE) %>% 
  as_tibble()

winprobs %<>% 
  mutate(pow_gy_bt = ifelse(is.na(pow_gy_bt),
                            0,
                            pow_gy_bt))


write_csv(winprobs, paste0(output, "log-abilities.csv"))



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

