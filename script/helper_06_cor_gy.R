library(tidyverse)
library(magrittr)
library(janitor)

list.files("data")

df <- "data/Et_diversity_panel_ALL-TRAITS.BLUPS.csv"
df %<>%
  read_csv() %>% 
  as_tibble(.name_repair = janitor::make_clean_names) %>% 
  rename(genotype = dna_code)



la <- "output/log-abilities/log-abilities.csv"
la %<>% 
  read_csv()

la %<>% 
  inner_join(. , df)

cor(la$pow_gy, la$gy_2012)
cor(la$pow_gy, la$gy_2013)

cor(la$pow_gy, rowMeans(la[,c("gy_2012","gy_2013")]))

cor(la$pow_rank, la$pow_gy)



la %>%
  select(genotype) %>%
  unique() %>%
  t() %>%
  as.vector() ->
  gene


 
df <- df[df$genotype %in% gene, ]




df <- df[c(3,which(grepl("gy", names(df))))]

df %<>% 
  select(-gy_2012, -gy_2013)


df %<>% 
  reshape2::melt(.) %>% 
  as_tibble()

df

df %<>% 
  mutate(id = as.integer(as.factor(variable))) %>% 
  arrange(id)


df %<>% 
  separate(., variable, c("a","year","region"), sep = "_") %>% 
  select(-a)





library(gosset)
library(PlackettLuce)


R <- to_rankings(df,
                 items = "genotype",
                 input = "value",
                 id = "id")


mod <- PlackettLuce(R)
summary(mod)
pow_gy_ep <- itempar(mod, log = TRUE)


la %<>% 
  mutate(pow_gy_ep = pow_gy_ep)


cor(la$gy_2013, pow_gy_ep)

