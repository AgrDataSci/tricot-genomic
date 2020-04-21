# Assess model performance and predictions trends
library("tidyverse")
library("magrittr")
library("janitor")
library("patchwork")
library("gosset")
library("PlackettLuce")

source("script/helper_00_functions.R")

sessioninfo::session_info()
# write session info
dir.create("script/session_info", recursive = TRUE, showWarnings = FALSE)
capture.output(sessioninfo::session_info(),
               file = "script/session_info/05_make_fig2.txt")

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
load("data/diversity.panel.data.gp.Rdata")
rm(snp.pos, geno, info, farm)

# climatology
# load("processing/climatology/climatology.rda")

# PL models
load("processing/plmodels/models.rda")
# load("processing/plmodels/models-100-fold.rda")

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
# .............................................
# .............................................
# Panel B, Reliability ####
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

# place predictions in a node
nodes <- vector()
for (i in seq_along(models)) {
  x <- predict(models[[i]], newdata = dt, type = "node")
  nodes <- cbind(nodes, x)
}

nodes <- ifelse(nodes[,1] == "2", "node 3",
                ifelse(nodes[,1] == "3" & nodes[,3] == "2", "node 4",
                       ifelse(nodes[,1] == "3" & nodes[,3] == "3", "node 5", NA)))

# get the reference variety from each plot
R <- dt$G
R <- R[1:length(R),,as.grouped_rankings = FALSE]

# define the number of best genotypes to compare with the reference
nbest <- 4

# run over obs
rel <- data.frame()
for(i in seq_len(n)) {
  
  # get the reference
  x <- R[i, ]
  x <- x[x!=0]
  f_impr <- improved[improved %in% names(x)]
  
  if (length(f_impr) == 0) next
  
  # get the bests but the reference 
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
                  reliability = as.numeric(re), 
                  node = rep(nodes[i], nbest))
  
  rel <- rbind(rel, d)
  
}

head(rel)

table(rel$ref)
table(rel$gen)

node_sel <- cbind(table(rel$gen[rel$node=="node 3"]),
                  table(rel$gen[rel$node=="node 4"]),
                  table(rel$gen[rel$node=="node 5"]))
dimnames(node_sel)

rel <- rel[rel$best < 4, ]
rel$best <- as.factor(ifelse(rel$best == 1, "1st",
                             ifelse(rel$best == 2,
                                    "2nd", "3rd")))

p2 <- ggplot(rel,
             aes(y = reliability, group = best, 
                 x = best, fill = best)) +
  geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values= c("#2166ac","#6baed6","#9ecae1")) + 
  labs(y = "Probability of outperforming",
       x = c(""),
       title = "B") +
  scale_y_continuous(limits = c(.8, .95),
                     breaks = seq(8, 10, 0.5)/10) +
  theme(axis.text.x = element_text(size=12, angle = 0,
                                   face="plain", colour = "black"),
        axis.text.y = element_text(size=12, angle = 0,
                                   hjust=1, vjust=0.5,
                                   face="plain", colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(3,3,1,3), "mm"),
        plot.title = element_text(size=16, 
                                  colour = "black", 
                                  face = "bold"),
        text = element_text(family = "sans"))

p2


# .....................................
# .....................................
# Panel A, Biplot ####
nodepca <- read.csv("output/reliability_yield_gain/node_gen.csv")
dt <- nodepca

nodesel <- data.frame(cbind(genotype = dimnames(node_sel)[[1]],
                 group = c(rep("3DB Cold tolerant", 3), rep("3DB Warm tolerant", 3))),
           stringsAsFactors = FALSE)

pca <- merge(nodepca, nodesel, by = "genotype", all.x = TRUE)
pca$group <- ifelse(pca$PC1 < - 100, "Currently recommended", pca$group) 
pca$group <- ifelse(is.na(pca$group), "Other genotypes", pca$group) 
pca$group <- factor(pca$group, levels = c("3DB Cold tolerant","3DB Warm tolerant",
                                          "Currently recommended","Other genotypes"))

p1<-
ggplot() +
  geom_point(data = pca, 
             aes(x = PC1, y = PC2, fill = group, 
                 color = group, shape = group), size = 1.5) +
  scale_shape_manual(values = c(15, 16, 17, 5)) +
  scale_color_manual(values= c("#91c43e","#9f49c3","#e0493d","gray40")) +
  labs(x = "PC1 (24.1%)",
       y = "PC2 (10.5%)",
       title = "A") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 0,
                                   face = "plain", colour = "black"),
        axis.text.y = element_text(size = 12, angle = 0,
                                   hjust = 1, vjust = 0.5,
                                   face = "plain", colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.25,0.9),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(3,3,1,3), "mm"),
        plot.title = element_text(size = 16, 
                                  colour = "black", 
                                  face = "bold"),
        text = element_text(family = "sans"))

p1
# .................................
# Panel C, Bar plot over the time series ####
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
# save(gen_gy, loc_gy, env_gy, baseline, gy, output,
#      file = paste0("output/reliability_yield_gain/", "predictions.rda"))
# load(paste0("output/reliability_yield_gain/", "predictions-100-fold.rda"))
# 
# year <- rep(1:15, each = length(loc_gy)*3)
# 
# year <- factor(year, levels = 1:15)

# gen_gy$year <- year
# env_gy$year <- year
# 
# yield <- bind_rows(env_gy, gen_gy)
# 
# yield$model <- factor(rep(c("ME","M3DB"), each = nrow(env_gy)), 
#                       levels = c("ME","M3DB"))
# 
# 
# yield %<>% 
#   mutate(gain = y / baseline -1)
# 
# e <- split(yield, yield$model)
# e <- lapply(e, function(x){
#   g <- glm(gain ~ year, data = x, family = poisson())
#   g <- sqrt(diag(vcov(g)))
#   as.vector(g)
# })
# 
# yield %<>% 
#   group_by(model, year) %>% 
#   summarise(gain = mean(gain)) %>% 
#   ungroup() %>% 
#   mutate(se = as.vector(unlist(e))) %>% 
#   filter(model == "M3DB")
# 
# yield %<>% 
#   mutate(se_min = gain - se,
#          se_max = gain + se)
# 
# write.csv(yield, "output/reliability_yield_gain/yield_gain.csv", row.names = FALSE)

yield <- read_csv("output/reliability_yield_gain/yield_gain.csv")

every_n_labeler <- function(n = 3) {
  function (x) {
    ind = ((1:length(x)) - 1) %% n == 0
    x[!ind] = ""
    return(x)
  }
}

p3 <- 
  ggplot(yield, aes(y = gain, x = year, fill = model)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = se_min, ymax = se_max), 
                width = 0.1,
                position = position_dodge(1)) +
  scale_x_discrete(breaks = seq(1,15,2)) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 3.5, 0.5)/10,
                     expand = c(0, 0)) +
  labs(y = "Increase in yield (%)",
       x = "Year",
       title = "C") +
  scale_fill_manual(values = c("#3182bd"), 
                      name = "") +
  theme(axis.text.x = element_text(size = 12, angle = 0,
                                   face="plain", colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, angle = 0,
                                   hjust=1, vjust=0.5,
                                   face="plain", colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size=12, colour="black"),
        legend.key = element_rect(colour = NA, fill = NA ,
                                  size=0.5,linetype = 1),
        legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linetype = "solid",
                                    fill = NA),
        plot.margin = unit(c(2,3,3,3), "mm"),
        plot.title = element_text(size=16, 
                                  colour = "black", 
                                  face = "bold"),
        text = element_text(family = "sans"))

p3

p <- p1 + (p2 / p3) + plot_layout(ncol = 2, widths = c(2, 1))

p

# save as svg
ggsave(paste0("manuscript/display_items/", "Fig2.svg"),
       plot = p,
       width = 20,
       height = 15,
       units = "cm")

# save as png
ggsave(paste0(paste0("manuscript/display_items/", "Fig2.png")),
       plot = p,
       width = 20,
       height = 15,
       dpi = 500,
       units = "cm")

