library("tidyverse")
library("magrittr")
library("gosset")
library("BradleyTerry2")



df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

output <- "output/exploring/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)


pow <- paste0(output, "probability_of_winning.csv")
pow %<>% 
  read_csv()

#.....................................
#.....................................
# farmer rank ####

# create grouped rankings 
R <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id")

B <- rank_binomial(R, drop.null = TRUE)

mod <- BTm(cbind(win1, win2), player1, player2, ~ genotype,
           id = "genotype", data = B, ref = "435ET_D")


summary(mod)

x <- coef(mod)

names(x) <- gsub("genotype","",names(x))

x <- tibble(genotype = names(x),
            pow_rank_bt = x)

pow %<>% 
  merge(. , x, by = "genotype", all.x = TRUE) %>% 
  as_tibble()


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


BY <- rank_binomial(YR, drop.null = TRUE)


mod <- BTm(cbind(win1, win2), player1, player2, ~ genotype,
           id = "genotype", data = BY)



x <- coef(mod)
names(x) <- gsub("genotype","",names(x))

x <- tibble(genotype = names(x),
            pow_gy_bt = x)


pow %<>% 
  merge(. , x, by = "genotype", all.x = TRUE) %>% 
  as_tibble()

pow %<>% 
  mutate(pow_gy_bt = ifelse(is.na(pow_gy_bt), 
                              0,
                              pow_gy_bt))

hist(pow$pow_gy_bt)

cor(pow$pow_gy, pow$pow_gy_bt)



write_csv(pow, paste0(output, "probability_of_winning.csv"))

