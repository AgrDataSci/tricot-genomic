# Assess model performance and predictions trends
library("tidyverse")
library("magrittr")
library("janitor")
library("patchwork")
library("gosset")
library("PlackettLuce")

source("script/helper_00_functions.R")

# ......................................
# .....................................
# Output dir
#output <- "output/reliability_yield_gain/"
output <- "manuscript/display_items/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

#.......................................
#.......................................
# Read data ####
df <- read_csv("data/durumwheat.csv")

# yield from station data
load("data/diversity.panel.data.gp.rda")
rm(snp.pos, geno, info, farm)

# climatology
load("processing/climatology/climatology.rda")

# PL models
load("processing/plmodels/models-100-fold.rda")

# varieties' names
items <- unique(df$genotype)
itemsG <- items[grepl("_D", items)]

# names of improved varieties
improved <- c("asassa", "ude", "hitosa")
# genotypes of improved varieties
improved <- unique(df[df$accession %in% improved, "genotype"])[[1]]
accession <- itemsG[!itemsG %in% improved]

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


# #.......................................
# #.......................................
# Plot results ####

# ..............................
# Pseudo R-squared boxplot ####
load("processing/plmodels/models.rda")
pr2 <- c(null$raw$estimators$MaxLik,
         gen$raw$estimators$MaxLik)

model <- factor(rep(c("MCB","M3DB"), 
                    each = length(null$raw$models)),
                levels = c("MCB","M3DB"))

model

pr2 <- tibble(value = pr2,
              model = model)

pr2

p1 <- ggplot(pr2,
             aes(y = value, group = model, x = model, fill = model)) +
  geom_boxplot(outlier.size = 0.5, size = 0.3, show.legend = FALSE) + 
  scale_fill_manual(values= c("#d73027","#2166ac")) +
  labs(y = bquote('Pseudo-R' ^2*''),
       x = "",
       title = "A") +
  scale_y_continuous(limits = c(0.4, 0.6),
                     breaks =  seq(40,60, 5)/100) + 
  scale_x_discrete(labels = c(MCB = "CB",
                              M3DB = "3DB")) +
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
p1

# .............................................
# .............................................
# Reliability plot ####

models <- gen$raw$models
dt <- gen$raw$data
n <- nrow(dt)

# predictions from cv models
pr <- array(NA, dim = c(n, length(items), 3), dimnames = list(1:n, items, 1:3))
for (i in seq_along(models)) {
  x <- predict(models[[i]], newdata = dt)
  pr[,,i] <- x
}
# average them 
pr <- apply(pr, c(1,2), mean)
pr <- pr[,c(improved, accession)]
pr <- t(apply(pr, 1, function(x){
  s <- sum(x)
  x <- x / s
}))

# get the reference variety from each plot
R <- dt$G
R <- R[1:length(R),,as.grouped_rankings = FALSE]

# define the number of best genotypes to compare with the reference
nbest <- 3

# run over obs
rel <- data.frame()
for(i in seq_len(n)) {
  
  # get the reference
  x <- R[i, ]
  x <- x[x!=0]
  f_impr <- improved[improved %in% names(x)]
  
  if (length(f_impr)==0) next
  
  # get the best three but the reference 
  y <- pr[i, ]
  p_impr <- y[f_impr]
  best <- y[!grepl(paste0(improved, collapse = "|"), names(y))]
  index_best <- sort(rank(best * -1))
  best <- best[names(index_best)[1:nbest]]
  
  re <- vector()
  for(j in seq_len(nbest)){
    r <- best[j] / (p_impr + best[j])
    re <- c(re, r)
  }
  
  d <- data.frame(id = rep(i, nbest),
                  best = 1:nbest,
                  ref = rep(names(p_impr), nbest),
                  gen = names(re),
                  reliability = as.numeric(re))
  
  rel <- rbind(rel, d)
  
}

head(rel)

table(rel$ref)
table(rel$gen)

rel$best <- as.factor(ifelse(rel$best == 1, "1st", ifelse(rel$best == 2, "2nd", "3rd")))

p2 <- ggplot(rel,
             aes(y = reliability, group = best, 
                 x = best, fill = best)) +
  geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values= c("#2166ac","#6baed6","#9ecae1")) + 
  labs(y = "Reliability",
       x = c(""),
       title = "B") +
  scale_y_continuous(limits = c(.8, .95),
                     breaks = seq(8, 10, 0.5)/10) +
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

p2

# .................................
# Bar plot over the time series ####
# # Gain 
# # gy by taking the best three varieties from Plackett-Luce models
# gen_gy <- NULL
# for(i in seq_along(climatology)){
#   print(i)
#   y <- predict_cv(gen, newdata = climatology[[i]])
#   y <- reliability_gy(y, gy, weights =  c(.7, .15, .15))
#   y <- data.frame(y = y, sample = paste0("S",i), r = 1:length(y))
#   gen_gy <- rbind(gen_gy, y)
# }
# 
# 
# 
# # 
# save(gen_gy, loc_gy, env_gy, baseline, gy, output,
#      file = paste0("output/reliability_yield_gain/", "predictions.rda"))
load(paste0("output/reliability_yield_gain/", "predictions-100-fold.rda"))

year <- rep(1:15, each = length(loc_gy)*3)

year <- factor(year, levels = 1:15)


gen_gy$year <- year
env_gy$year <- year

yield <- bind_rows(env_gy, gen_gy)

yield$model <- factor(rep(c("ME","M3DB"), each = nrow(env_gy)), 
                      levels = c("ME","M3DB"))


yield %<>% 
  mutate(gain = y / baseline -1)

e <- split(yield, yield$model)
e <- lapply(e, function(x){
  g <- glm(gain ~ year, data = x, family = poisson())
  g <- sqrt(diag(vcov(g)))
  as.vector(g)
})

yield %<>% 
  group_by(model, year) %>% 
  summarise(gain = mean(gain)) %>% 
  ungroup() %>% 
  mutate(se = as.vector(unlist(e))) %>% 
  filter(model == "M3DB")

yield %<>% 
  mutate(se_min = gain - se,
         se_max = gain + se)

p3 <- 
  ggplot(yield, aes(y = gain, x = year, fill = model)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = se_min, ymax = se_max), 
                width = .1,
                position = position_dodge(1)) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 3.5, 0.5)/10 ) +
  labs(y = "Increase in yield (%)",
       x = "Year",
       title = "C") +
  scale_fill_manual(values = c("#3182bd"), 
                      name = "") +
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
        #legend.background = element_rect(colour = NA, fill = NA),
        #legend.justification = c(1, 0), 
        legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(2,3,3,3), "mm"),
        plot.title = element_text(size=16, 
                                  colour = "black", 
                                  face = "bold"))

p3

p <- (p1 | p2) / p3

p

# save as svg
ggsave(paste0("manuscript/display_items/", "Fig2.svg"),
       plot = p,
       width = 15,
       height = 15,
       units = "cm")

# save as png
ggsave(paste0(paste0("manuscript/display_items/", "Fig2.png")),
       plot = p,
       width = 15,
       height = 15,
       units = "cm")
