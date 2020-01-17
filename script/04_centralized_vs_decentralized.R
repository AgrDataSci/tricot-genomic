# Compare centralized vs decentralized

library("tidyverse")
library("magrittr")
library("janitor")
library("PlackettLuce")
library("gosset")

#.................................................
#.................................................
# Data ####
load("data/diversity.panel.data.gp.Rdata")

# farmer rankings
df <- read_csv("data/durumwheat.csv")
rm(snp.pos, geno, info)

#.....................................
#.....................................
#.....................................
# Overall appreciation ####

# get oa from decentralized trials
items <- unique(df$genotype)
items <- items[grepl("_D", items)]

# keep only genotyped items
keep <- df$genotype %in% items

oa <- 
  df %>% 
  filter(keep)

# remove rankings with less than 2 items
keep <-
  oa %>%
  group_by(id) %>%
  summarise(keep = n_distinct(genotype) >= 2) %>%
  filter(keep)

keep <- oa$id %in% keep$id

oa %<>% 
  filter(keep)


# .....................................
# 'farm' has the data of ao for centralized trials 
keep <- which(grepl("OVERALL_", names(farm)) & !grepl("_L", names(farm)))
keep <- c(1:8, keep)

oa2 <- farm[, keep]

oa2 %<>% 
  as_tibble(.name_repair = make_clean_names) %>% 
  rename(genotype = id)

# keep only genotyped items
keep <- oa2$genotype %in% items

oa2 %<>% 
  filter(keep)

oa2 %<>% 
  select(-plot, -blk, -row, -col)


# rename columns 
target <- which(grepl("_overall_", names(oa2)))
names(oa2)[target] <- paste0("oa", seq_along(target))

oa2 %<>% 
  pivot_longer(cols = starts_with("oa"),
               names_to = "farmer", 
               names_prefix = "oa",
               values_to = "oa",
               values_drop_na = TRUE)


# create an id
oa2 %<>% 
  mutate(id = paste0(year, locality, rep, farmer)) %>% 
  mutate(id = as.integer(as.factor(id)))





oa <- rank_numeric(data = oa,
                   items = "genotype",
                   input = "farmer_rank",
                   id = "id", 
                   ascending = TRUE)

oa <- oa[1:length(oa), , as.rankings = FALSE]


G2 <- rank_numeric(oa2, 
                   "genotype", 
                   "oa", 
                   "id")

G2 <- G2[1:length(G2), , as.rankings = FALSE]


n <- dim(oa)[1]
kendall <- rep(NA, n)

n2 <- dim(G2)[[1]]

for (j in seq_len(n)) {
  
  f <- oa[j, oa[j, ] != 0]
  
  k <- rep(NA, n2)
  
  for (i in seq_len(n2)) {
  
    g2_i <- G2[i, names(f)]
    
    if (sd(g2_i) == 0) next
      
    k_i <- try(
      kendallTau(f,
                 g2_i)$kendallTau
      , silent = TRUE)
    
    k[i] <- k_i
    
  }
  
  k <- mean(k, na.rm = TRUE)
  
  kendall[j] <- k

}

mean(kendall)

#.....................................
#.....................................
#.....................................
# GRAIN YIELD ####

#...........................
#...........................
# grain yield from decentralized
# remove NAs in grain yield
df %>% 
  filter(!is.na(gy_gm)) -> 
  gy


gy %<>% 
  group_by(id) %>% 
  distinct(gy_gm, .keep_all = TRUE) %>% 
  ungroup()

keep <- gy$genotype %in% items

gy %<>% 
  filter(keep)


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
gy <- gy[id, ]

gy

gy <- rank_numeric(gy,
                   "genotype",
                   "gy_gm",
                   "id")

gy <- gy[1:length(gy),,as.rankings = FALSE]


# met has the agronomic data from centralized trials
head(met)

met %>% 
  select(LOCATION, YEAR, REP, ID, GY) %>% 
  as_tibble(.name_repair = make_clean_names) %>% 
  rename(genotype = id) ->
  gy2
  

gy2 %<>% 
  mutate(id = as.integer(as.factor(paste0(location, year, rep))))


keep <- gy2$genotype %in% items


gy2 %<>% 
  filter(keep)

gy2 <- rank_numeric(gy2,
                    "genotype",
                    "gy", 
                    "id")

gy2 <- gy2[1:length(gy2),, as.rankings = FALSE]



n <- dim(gy)[1]
kendall <- rep(NA, n)

n2 <- dim(gy2)[[1]]

for (j in seq_len(n)) {
  
  f <- gy[j, gy[j, ] != 0]
  
  k <- rep(NA, n2)
  
  for (i in seq_len(n2)) {
    
    g2_i <- gy2[i, names(f)]
    
    if (sd(g2_i) == 0) next
    
    k_i <- try(
      kendallTau(f,
                 g2_i)$kendallTau
      , silent = TRUE)
    
    k[i] <- k_i
    
  }
  
  k <- mean(k, na.rm = TRUE)
  
  kendall[j] <- k
  
}

mean(kendall)

boxplot(kendall)

