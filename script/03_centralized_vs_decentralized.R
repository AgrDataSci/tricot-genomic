# Compare centralized vs decentralized
# Here we compare the rankings provided by the centralized trials
# with the rankings provided by centralized using kendall correlation

# Packages
library("tidyverse")
library("magrittr")
library("janitor")
library("PlackettLuce")
library("gosset")
library("ggplot2")

#.................................................
#.................................................
# Data ####
load("data/diversity.panel.data.gp.Rdata")
rm(snp.pos, geno, info)

# farmer rankings
df <- read_csv("data/durumwheat.csv")

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

oa2 %<>% 
  group_by(locality, genotype, year) %>% 
  summarise(oa = mean(oa)) %>% 
  ungroup()


oa2 %<>% 
  mutate(id = paste0(year, locality)) %>% 
  mutate(id = as.integer(as.factor(id)))
  


g1 <- rank_numeric(data = oa,
                   items = "genotype",
                   input = "farmer_rank",
                   id = "id")

g1 <- g1[1:length(g1), , as.rankings = FALSE]


g2 <- rank_numeric(oa2,
                   items = "genotype", 
                   input = "oa", 
                   id = "id")

g2 <- g2[1:length(g2), , as.rankings = FALSE]

n <- dim(g1)[1]
kendall <- rep(NA, n)

n2 <- dim(g2)[[1]]

for (j in seq_len(n)) {
  
  # keep only the values with rankings
  g1_i <- g1[j, g1[j, ] != 0]
  
  # create a vector to keep all kendall taus
  k <- rep(NA, n2)
  
  # run over centralised trails
  for (i in seq_len(n2)) {
  
    # keep only the same genotypes as found in g1_i
    g2_i <- g2[i, names(g1_i)]
    
    # if all ties then jump to the next
    if (sd(g2_i) == 0) next
    
    # compute kendall tau  
    k_i <- try(
      kendallTau(g1_i,
                 g2_i)$kendallTau
      , silent = TRUE)
    
    k[i] <- k_i
    
  }
  
  k <- mean(k, na.rm = TRUE)
  
  kendall[j] <- k

}

mean(kendall)

boxplot(kendall)

oa_kendall <- kendall

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

gy <- rank_numeric(gy,
                   "genotype",
                   "gy_gm",
                   "id")

gy <- gy[1:length(gy), ,as.rankings = FALSE]


# met has the agronomic data from centralized trials
head(met)

met %>% 
  select(LOCATION, YEAR, REP, ID, GY) %>% 
  as_tibble(.name_repair = make_clean_names) %>% 
  rename(genotype = id) ->
  gy2
  

gy2 %<>% 
  group_by(location, year, genotype) %>% 
  summarise(gy = mean(gy)) %>% 
  ungroup()

gy2 %<>% 
  mutate(id = as.integer(as.factor(paste0(location, year))))


gy2 <- rank_numeric(gy2,
                    "genotype",
                    "gy", 
                    "id")

gy2 <- gy2[1:length(gy2),, as.rankings = FALSE]


n <- dim(gy)[1]
kendall <- rep(NA, n)

n2 <- dim(gy2)[[1]]

for (j in seq_len(n)) {
  
  gy_i <- gy[j, gy[j, ] != 0]
  
  k <- rep(NA, n2)
  
  for (i in seq_len(n2)) {
    
    gy2_i <- gy2[i, names(gy_i)]
    
    if (sd(gy2_i) == 0) next
    
    k_i <- try(
      kendallTau(gy_i,
                 gy2_i)$kendallTau
      , silent = TRUE)
    
    k[i] <- k_i
    
  }
  
  k <- mean(k, na.rm = TRUE)
  
  kendall[j] <- k
  
}

mean(kendall)

boxplot(kendall)

gy_kendall <- kendall

round(mean(oa_kendall), 4)
round(mean(gy_kendall), 4)

kendall <- data.frame(model = c(rep("OA", length(oa_kendall)), 
                                rep("GY", length(gy_kendall))),
                      kendall = c(oa_kendall,
                                  gy_kendall),
                      stringsAsFactors = FALSE)

head(kendall)

p <-
ggplot(kendall,
       aes(y = kendall, group = model, x = model)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(y = "",
       x = "",
       title = "") +
  theme(axis.text.x = element_text(size=14, angle = 0,
                                   face="plain", colour = "black"),
        axis.text.y = element_text(size=14, angle = 0,
                                   hjust=1, vjust=0.5,
                                   face="plain", colour = "black"),
        axis.title.y = element_text(size=14, colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(3,3,1,3), "mm"),
        plot.title = element_text(size=16, 
                                  colour = "black", 
                                  face = "bold"))

# save as svg
ggsave("output/SI/kendall_cor_cent_vs_decent.svg",
       plot = p,
       width = 15,
       height = 15,
       units = "cm")

# save as png
ggsave("output/SI/kendall_cor_cent_vs_decent.png",
       plot = p,
       width = 15,
       height = 15,
       units = "cm")
# 

