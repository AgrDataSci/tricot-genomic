# Bradley-Terry model with player specific abilities
# Kaue de Sousa
# INN University / Bioversity International
# Updated 12Jun2019

# Load required packages 
library("tidyverse")
library("magrittr")
library("gosset")
library("BradleyTerry2")
library("caret")



pairwise_contests <- function(object) {
  
  if(class(object) == "rankings") {
    object <- to_paircomp(object)
  }
  
  if(class(object) == "grouped_rankings") {
    object <- to_paircomp(object)
  }
  
  if(class(object) != "paircomp") {
    stop("Unknown class for object")
  }
  
  P <- object
  
  P <- as.matrix(P)
  
  contests <- apply(P, 1, function(X) {
    
    X <- X[!is.na(X)]
    
    players <- as.matrix(do.call(rbind, strsplit(names(X), ":")))
    
    X <- as.vector(ifelse(X == -1, 2, 1))
    
    w <- cbind(1:length(X), X)
    
    winner <- players[w]
    
    X <- ifelse(X == 2, 1, 2)
    
    l <- cbind(1:length(X), X)
    
    loser <- players[l]
    
    return(cbind(winner, loser))
  })
  
  
  contests <- do.call(rbind, contests)
  
  if (any(contests[,1] == contests[,2])) {
    warning("winner and loser has equal values")
  }
  
  
  contests <- data.frame(contests, stringsAsFactors = FALSE)
  
  items <- unique(unlist(contests[,1:2]))
  
  contests$winner <- factor(contests$winner, levels = items)
  contests$loser <- factor(contests$loser, levels = items)
  
  return(contests)
  
}

# create the output path for this analysis
output <- "output/bradley-terry_gene/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# read the data
df <- "data/durumwheat.csv"
df %<>% 
  read_csv()

#.....................................
#.....................................
# farmer rank ####

# create the contests object 
R <- to_rankings(data = df,
                 items = "genotype",
                 input = "farmer_rank",
                 id = "id")


# get the dataframe of contests between items
contests <- pairwise_contests(R)

items <- sort(as.character(unique(unlist(contests))))

df %>% 
  distinct(accession, .keep_all = TRUE) %>% 
  select(accession, genotype) ->
  labels

labels <- as.vector(t(labels[match(items, labels$genotype), 1]))

#...........................................
#...........................................
# create a matrix with predictors ####

# read genotypic data
load("data/genotypic.data.durum.wheat.rda")

gene <- imputedcl
rm(imputedcl)
dimnames(gene)

# keep only the items used in this analysis
gene <- gene[dimnames(gene)[[1]] %in% items, ]

# remove those with no variance
drop <- caret::nearZeroVar(gene)

gene <- gene[,-drop]

pca <- prcomp(gene, center = TRUE,scale. = TRUE)

predictors <- pca$x[,c(1:4)] 

predictors %<>% 
  as.data.frame()


names(predictors) <- paste0("gene_", names(predictors))


#..........................................
#..........................................
# add genotype traits from field plots ####
load("data/blup.diversity.panel.rda")

drop <- rownames(phenos) %in% items

phenos %<>%
  select(GY, DF, TGW, SPS, SPL)

phenos <- phenos[drop, ]

predictors <- cbind(predictors, phenos)

# combine datasets with contests and explanatory variables
durum <- list(contests = contests, 
              predictors = predictors)


#..........................................
#..........................................
# fit the model ####

# model without explanatory variables
result <- rep(1, nrow(contests))
modn <- BTm(result, winner, loser, 
            family = binomial(link = "probit"),
            data = durum$contests) 


summary(modn)



# model with variables
mod <- BTm(result, winner, loser,
           ~ gene_PC1[..] + GY[..] + DF[..] + SPL[..] +
             SPS[..] + TGW[..],
           data = durum)

summary(mod)

AIC(modn)

AIC(mod)


plot(BTabilities(modn)[,1], las = 2)
axis(1, at=1:length(labels) , labels=labels, las = 2)

