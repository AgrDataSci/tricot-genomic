# Model scenarios with explanatory variables 
# to predict variety performance using PlackettLuce models

library("tidyverse")
library("magrittr")
library("PlackettLuce")
library("gosset")
library("patchwork")
library("caret")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/04_PL_model_session_info.txt")

#.................................................
#.................................................
# Data ####
list.files("data")

# farmer rankings
df <- read_csv("data/durumwheat.csv")

# environmental indices
ind <- read_csv("data/environmental_indices.csv")

# genotypic data
load("data/genotypic.data.durum.wheat.rda")
gene <- imputedcl
rm(imputedcl, hmp)

#.................................................
#.................................................
# PlackettLuce rankings ####

G <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id",
                  group = TRUE)


R <- rank_numeric(data = df,
                  items = "genotype",
                  input = "farmer_rank",
                  id = "id",
                  group = FALSE)


#.................................................
#.................................................
# normal prior coefficients ####
mod <- PlackettLuce(G)

items <- names(coef(mod))

genitems <- items[items %in% dimnames(gene)[[1]]]

keep <- dimnames(gene)[[1]] %in% genitems

# keep only the items used in this analysis
gene <- gene[keep, ]

# remove those with no variance
drop <- nearZeroVar(gene)
gene <- gene[, -drop]

set.seed(765)
am <- matrix(runif(length(items) * ncol(gene), -1e-1, 1e-1), 
             nrow = length(items), 
             ncol = dim(gene)[[2]],
             dimnames = list(items,
                             dimnames(gene)[[2]]))

am[genitems,] <- gene

am <- t(am)

am <- cov(am)

chol.default(am)

prior <- list(mu = as.vector(coef(mod)),
              Sigma = am)

#.................................................
#.................................................
# Forward selection ####
output <- "processing/plmodels/"
dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

# cross validation parameters
n <- length(G)
minsize <- round(n * 0.30, -2)
npseudo <- 5
alpha <- 0.01
bonferroni <- FALSE
gamma <- TRUE
folds <- as.numeric(as.factor(ind$year))
k <- max(folds)


# select variables
keep <- which(grepl("DT_|NT_|MLDS_rep|lon|lat|xy|yx", names(ind)))

ind <- ind[, keep]

# put all together as mydata
mydata <- cbind(G, ind)

#.................................................
#.................................................
# Model with OA ####
gen <- crossvalidation(G ~ minNT_veg + maxNT_rep,
                       data = mydata,
                       k = k,
                       folds = folds,
                       minsize = minsize,
                       alpha = alpha,
                       normal = prior,
                       gamma = gamma)


# Now using GY
#...........................
# grain yield from decentralized
# remove NAs in grain yield
gy <- rank_numeric(df,
                   "genotype",
                   "gy_gm",
                   "id")

keep <- !is.na(gy)

gy <- gy[keep, ]

gy_ind <- ind[keep, ]

GY <- group(gy, index = 1:length(gy))

dataGY <- cbind(GY, gy_ind)

folds_gy <- folds[keep]

gen_gy <- crossvalidation(GY ~ minNT_veg + maxNT_rep,
                          data = dataGY,
                          k = k,
                          folds = folds_gy,
                          minsize = minsize,
                          alpha = alpha,
                          normal = prior)



save(gen, gen_gy,
     file = paste0(output, "models.rda"))

#...................................
#...................................
# Plot tree ####
plt <- pltree(formula(gen$raw$call),
              data = mydata,
              alpha = alpha,
              minsize = 100,
              normal = prior,
              gamma = TRUE)


# Extract ids from terminal nodes
node_id <- partykit::nodeids(plt, terminal = TRUE)

# get node information
nodes <- list()
for (i in seq_along(node_id)) {
  nodes[[i]] <- plt[[ node_id[i] ]]$node$info$object
}

# get number of observers in each node
nobs <- list()
for (i in seq_along(node_id)) {
  nobs[[i]] <- as.integer(plt[[ node_id[i] ]]$node$info$nobs) 
}

# get item parameters from model
coeffs <- lapply(nodes, psychotools::itempar)

# get estimates from item parameters using qvcalc
coeffs <- lapply(coeffs, qvcalc::qvcalc)

# extract dataframes with estimates
coeffs <- lapply(coeffs, function(X){
  df <- X[]$qvframe }
)

# get item names
items <- rownames(coeffs[[1]])
items <- items[grepl("_D", items)]
# Add limits in error bars and item names
coeffs <- lapply(coeffs, function(X){
  X <- X[items, ]
  s <- sum(X$estimate)
  X$estimate <- X$estimate / s
  X <- within(X, {
    bmin <- X$estimate-(X$quasiSE)
    bmax <- X$estimate+(X$quasiSE)
    items <- items
  })
  
  X$bmax <- ifelse(X$bmax > 1, 0.991, X$bmax)
  
  X$bmin <- ifelse(X$bmin < 0, 0.001, X$bmin)
  return(X)
})

# Add node information and number of observations
for (i in seq_along(node_id)) {
  coeffs[[i]] <- within(coeffs[[i]], {
    nobs <- nobs[[i]]
    node <- node_id[i]}
  )
}

# Get max and min values for the x axis in the plot
xmax <- round(max(do.call(rbind, coeffs)$bmax, na.rm = TRUE) + 0.01, digits = 4)
xmin <- round(min(do.call(rbind, coeffs)$bmin, na.rm = TRUE), digits = 4)

# Check font size for axis X and Y, and plot title
s.title <- 13
s.axis <- 15
labels <- coeffs[[1]]$items

# Plot winning probabilities
plots <- list()
for(i in seq_along(coeffs)) {
  
  X <- coeffs[[i]]
  
  p <- ggplot2::ggplot(X, ggplot2::aes(x = estimate, y = labels)) +
    ggplot2::geom_vline(xintercept = 1/length(items), 
                        colour = "#E5E7E9", size = 0.8) +
    ggplot2::geom_point(pch = 21, size = 2, 
                        fill = "black",colour = "black") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = bmin,
                                         xmax = bmax),
                            colour="black", height = 0.2) +
    ggplot2::scale_x_continuous(limits = c(0, xmax)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL ,
                  title = paste0("Node ", X$node[1], " (n= ", X$nobs[1], ")")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = s.title),
                   axis.text.x = ggplot2::element_text(size = s.axis, angle = 0,
                                                       hjust = 0.5, vjust = 1, face = "plain",
                                                       colour = "black"),
                   axis.text.y = ggplot2::element_text(size = s.axis, angle = 0,
                                                       hjust = 1, vjust = 0.5, face = "plain",
                                                       colour = "black"),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(colour = "black", size = 1),
                   axis.ticks = ggplot2::element_line(colour = "black", size = 0.5),
                   axis.ticks.length = grid::unit(0.3, "cm"))
  
  plots[[i]] <- p
  
}

plots

p <- plots[[1]] | plots[[2]] | plots[[3]]
p
ggsave("output/SI/pltree.svg",
       plot = p,
       width = 25,
       height = 25,
       units = "cm")

svg("output/SI/plt.svg", width = 7, height = 8, pointsize = 15)
plot(plt)
dev.off()

# Plot R-squared
pr2 <- c(as.vector(unlist(gen_gy$raw$estimators[,7])),
         as.vector(unlist(gen$raw$estimators$MaxLik)))

model <- factor(rep(c("GY","OA"),
                    each = length(gen$raw$models)),
                levels = c("GY","OA"))

model
pr2

pr2 <- tibble(value = pr2,
              model = model)

pr2

p1 <- ggplot(pr2,
             aes(y = value, group = model, x = model, fill = model)) +
  geom_boxplot(outlier.size = 0.5, size = 0.3, show.legend = FALSE) +
  scale_fill_manual(values= c("#d73027","#2166ac")) +
  labs(y = bquote('Pseudo-R' ^2*''),
       x = "",
       title = "") +
  scale_y_continuous(limits = c(0.35, 0.6),
                     breaks =  seq(35,60, 5)/100) +
  theme(axis.text.x = element_text(size = 12, angle = 0,
                                   face = "plain", colour = "black"),
        axis.text.y = element_text(size = 12, angle = 0,
                                   hjust = 1, vjust = 0.5,
                                   face = "plain", colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
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

# save as svg
ggsave(paste0("manuscript/display_items/SI/", "FigS13.png"),
       plot = p1,
       width = 10,
       height = 10,
       units = "cm")

# #......................................
# #......................................
# # now predict OA >> GY and GY >> OA
# source("script/helper_00_functions.R")
# 
# pred_gen_gy <- predict_cv(gen_gy, newdata = mydata)
# 
# OA <- G[1:length(G), , as.grouped_rankings = FALSE]
# 
# kOA <- rep(NA, nrow(OA))
# 
# for(i in seq_along(kOA)){
#   kOA[i] <- kendallTau(pred_gen_gy[i, ], OA[i, ])[[1]]
# }
# 
# mean(kOA, na.rm = TRUE)
# sd(kOA, na.rm = TRUE)
# 
# pred_oa_gy <- predict_cv(gen, newdata = dataGY)
# 
# GY <- GY[1:length(GY), , as.grouped_rankings = FALSE]
# 
# kGY <- rep(NA, nrow(GY))
# 
# for(i in seq_along(kGY)){
#   kGY[i] <- kendallTau(pred_oa_gy[i, ], GY[i, ])[[1]]
# }
# 
# mean(kGY)
# sd(kGY)
