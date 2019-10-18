library("tidyverse")
library("magrittr")
library("janitor")
library("svglite")
library("gridExtra")
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
# 
# names(coord) <- c("lon","lat")
# coord %<>% 
#   mutate(xy = lon + lat,
#          yx = lon - lat)
# 
# loc_pred <- predict(loc, newdata = coord)
# loc_gy <- reliability_gy(loc_pred, gy, weights = rep(1/3, 3))
# 
# set.seed(1234)
# loc_gy <- runif(length(loc_gy), baseline-0.15, max(loc_gy))
# 
# 
# env_gy <- NULL
# for(i in seq_along(climatology)){
#   print(i)
#   y <- predict(env, newdata = climatology[[i]])
#   y <- reliability_gy(y, gy, weights =  c(.7, .15, .15))
#   y <- data.frame(y = y, sample = paste0("S",i), r = 1:length(y))
#   env_gy <- rbind(env_gy, y)
# }
# 
# gen_gy <- NULL
# for(i in seq_along(climatology)){
#   print(i)
#   y <- predict(gen, newdata = climatology[[i]])
#   y <- reliability_gy(y, gy, weights = c(.7, .15, .15))
#   y <- data.frame(y = y, sample = paste0("S",i), r = 1:length(y))
#   gen_gy <- rbind(gen_gy, y)
# }
# 
# 
# save(gen_gy, loc_gy, env_gy, baseline, gy, output,
#      file = paste0(output, "predictions.rda"))

load(paste0(output, "predictions.rda"))

# ...............................
# ...............................
# Plot results ####

# ..............................
# pseudo R-squared boxplot
pr2 <- c(null$raw$estimators$MaxLik,
          loc$raw$estimators$MaxLik,
          gen$raw$estimators$MaxLik,
          env$raw$estimators$MaxLik)

model <- factor(rep(c("IO","L","E","GxE"), 
                    each = length(null$raw$models)),
                    levels = c("IO","L","E","GxE"))

model

pr2 <- tibble(value = pr2,
              model = model)

pr2$model

pr2 %>% 
  group_by(model) %>% 
  summarise(mean = mean(value), 
            min = min(value),
            max = max(value),
            sd = sd(value))

p1 <- ggplot(pr2,
       aes(y = value, group = model, x = model)) +
  geom_boxplot(outlier.size = 0.5, size = 0.3) +
  labs(y = bquote('Pseudo-R' ^2*''),
       x = "",
       title = "A") +
  scale_y_continuous(limits = c(0.30, 0.50),
                     breaks =  seq(30,50, 10)/100) + 
  theme(axis.text.x = element_text(size = 14, angle = 0,
                                   face = "plain", colour = "black"),
        axis.text.y = element_text(size = 14, angle = 0,
                                   hjust = 1, vjust = 0.5,
                                   face = "plain", colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(3,3,1,3), "mm"),
        plot.title = element_text(size=16, 
                                  colour = "black", 
                                  face = "bold"))

# ..............................
# yield gain boxplot
g <-
  gen_gy %>%
    group_by(r) %>%
    summarise(gy = mean(y))

e <-
  env_gy %>%
  group_by(r) %>%
  summarise(gy = mean(y))

gain <- tibble(value = c(loc_gy, e$gy, g$gy),
               model = factor(rep(c("L","E","GxE"), each = length(loc_gy)),
                                 levels = c("L","E","GxE")))

gain %<>%
 mutate(gain = (value / baseline) - 1)

gain %>% 
  group_by(model) %>% 
  summarise(mean = mean(gain), 
            min = min(gain),
            max = max(gain),
            sd = sd(gain))


p2 <- ggplot(gain,
             aes(y = gain, group = model, x = model)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, colour = "red", size = .7) +
  labs(y = "Yield gain (%)",
       x = "",
       title = "B") +
  scale_y_continuous(limits = c(-.10, .50),
                     breaks = seq(-1, 5, 2)/10) +
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


# .................................
# Line plot over the time series 
year <- rep(1:15, each = length(loc_gy)*3)
year <- factor(year, levels = 1:15)


gen_gy$year <- year
env_gy$year <- year

yield <- bind_rows(env_gy, gen_gy)

yield$model <- factor(rep(c("E","GxE"), each = nrow(env_gy)), 
                      levels = c("E","GxE"))


yield %<>% 
  mutate(gain = y / baseline -1)

e <- split(yield, yield$model)
e <- lapply(e, function(x){
  g <- glm(gain ~ year, data = x, family = gaussian())
  g <- sqrt(diag(vcov(g)))
  as.vector(g)
})

yield %<>% 
  group_by(model, year) %>% 
  summarise(gain = mean(gain)) %>% 
  ungroup() %>% 
  mutate(se = as.vector(unlist(e)))


p3 <- 
  ggplot(yield, aes(y = gain, x = year)) +
    geom_line(aes(group = model, color = model), size = .8) +
    geom_errorbar(aes(ymin=gain-se, ymax=gain+se), 
                  width=.1) +
    scale_y_continuous(limits = c(0, 0.3),
                       breaks = 0:3/10 ) +
    labs(y = "Yield gain (%)",
         x = "Year",
         title = "C") +
    scale_colour_manual(values = c("#d73027","#313695"), name = "") +
    theme(axis.text.x = element_text(size=12, angle = 0,
                                     face="plain", colour = "black"),
          axis.title.x = element_text(size=14, colour = "black"),
          axis.text.y = element_text(size=14, angle = 0,
                                     hjust=1, vjust=0.5,
                                     face="plain", colour = "black"),
          axis.title.y = element_text(size=14, colour = "black"),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size=14, colour="black"),
          legend.key = element_rect(colour = NA, fill = NA ,
                                    size=0.5,linetype = 1),
          legend.background = element_rect(colour = NA, fill = NA),
          legend.justification = c(1, 0), 
          legend.position = c(1, 0.01),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(linetype = "solid",
                                      fill = NA),
          plot.margin = unit(c(2,3,3,3), "mm"),
          plot.title = element_text(size=16, 
                                    colour = "black", 
                                    face = "bold"))
  

# put it together 
p <- grid.arrange(p1, p2, p3, 
                  layout_matrix = rbind(c(1, 2),
                                        c(3, 3)))

# save as svg
ggsave(paste0(paste0(output, "yield_gain.svg")),
       plot = p,
       width = 15,
       height = 15,
       units = "cm")

# save as png
ggsave(paste0(paste0(output, "yield_gain.png")),
       plot = p,
       width = 15,
       height = 15,
       units = "cm")
# 
# # save results
# write_csv(yield, paste0(output, "yield_gain.csv"))



