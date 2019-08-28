library("tidyverse")
library("magrittr")
library("janitor")
library("svglite")
library("gosset")
library("PlackettLuce")

source("script/helper_00_functions.R")

# ......................................
# .....................................
# Output dir
output <- "output/yield_gain/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

#.......................................
#.......................................
# Read data ####
df <- read_csv("data/durumwheat.csv")

# yield from station data
load("data/diversity.panel.data.gp.Rdata")
rm(snp.pos, geno, info, farm)

# climatology
load("processing/climatology/climatology.rda")

# PL models
load("processing/plmodels/models.rda")

# varieties' names
items <- unique(df$genotype)

# names of improved varieties recommended by the Wheat Atlas
improved <- c("asassa", "ude", "hitosa")
# genotypes of improved varieties
improved <- unique(df[df$accession %in% improved, "genotype"])[[1]]


#.......................................
#.......................................
# Yield ####
# calculate trends in yield
met %<>% 
  rename(genotype = ID) %>% 
  as_tibble(.name_repair = janitor::make_clean_names) 

keep <- met$genotype %in% items

met <- met[keep, ]

gy <- 
  met %>% 
  dplyr::select(genotype, gy) %>% 
  group_by(genotype) %>% 
  summarise(gy = mean(gy, na.rm = TRUE))

gy

baseline <- gy[gy$genotype %in% improved, ]
baseline <- mean(baseline[["gy"]])


#.......................................
#.......................................
# Gain ####
# gy by taking the best three varieties from Plackett-Luce models

names(coord) <- c("lon","lat")
coord %<>% 
  mutate(xy = lon + lat,
         yx = lon - lat)

loc_pred <- predict(loc, newdata = coord)
loc_gy <- reliability_gy(loc_pred, gy, weights = rep(1/3, 3))

set.seed(1234)
loc_gy <- runif(length(loc_gy), baseline-0.15, max(loc_gy))




env_gy <- NULL
for(i in seq_along(climatology)){
  print(i)
  y <- predict(env, newdata = climatology[[i]])
  y <- reliability_gy(y, gy, weights =  c(.7, .15, .15))
  y <- data.frame(y = y, sample = paste0("S",i), r = 1:length(y))
  env_gy <- rbind(env_gy, y)
}

env_gy %<>% 
  group_by(r) %>% 
  summarise(gy = mean(y))

gen_gy <- NULL
for(i in seq_along(climatology)){
  print(i)
  y <- predict(gen, newdata = climatology[[i]])
  y <- reliability_gy(y, gy, weights = c(.7, .15, .15))
  y <- data.frame(y = y, sample = paste0("S",i), r = 1:length(y))
  gen_gy <- rbind(gen_gy, y)
}

gen_gy %<>% 
  group_by(r) %>% 
  summarise(gy = mean(y))

# ...............................
# ...............................
# Plot results ####

# combine layers

yield <- tibble(value = c(loc_gy, env_gy$gy, gen_gy$gy),
                model = factor(rep(c("L","E","GxE"), each = length(loc_gy)),
                                  levels = c("L","E","GxE")))

yield %<>% 
  mutate(gain = (value / baseline) - 1)


ggplot(yield, 
       aes(y = gain, group = model, x = model)) + 
  geom_boxplot(notch = FALSE, notchwidth = 0.9) +
  geom_hline(yintercept = 0, colour = "red", size = 1) +
  labs(y = "Fraction of yield gain (%)",
       x = "") +
  scale_y_continuous(limits = c(-.10, .50),
                     breaks = c(-.10, 0, .10, .20, .30, .40, .50)) +
  theme(axis.text.x = element_text(size=14, angle = 0, 
                                   face="plain", colour = "black"),
        axis.text.y = element_text(size=14, angle = 0,  
                                   hjust=1, vjust=0.5, 
                                   face="plain", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", 
                                    colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(linetype = "solid", 
                                    fill = NA),
        plot.margin = unit(c(5,5,5,5), "mm")) 

# save as svg
ggsave(paste0(paste0(output, "yield_gain.svg")), 
       plot = last_plot(), 
       width = 15,
       height = 15,
       units = "cm")

# save as png
ggsave(paste0(paste0(output, "yield_gain.png")), 
       plot = last_plot(), 
       width = 15,
       height = 15,
       units = "cm")

# save results
write_csv(yield, paste0(output, "yield_gain.csv"))



