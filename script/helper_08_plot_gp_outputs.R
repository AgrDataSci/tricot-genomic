library("tidyverse")
library("magrittr")
library("reshape2")
library("svglite")
library("gridExtra")

# read data produced in script 1.rrblup.R
load("data/GP.CSdata.Rdata")

#['#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac']
#['#762a83','#af8dc3','#e7d4e8','#d9f0d3','#7fbf7b','#1b7837']

# produce some plotting ####
gptraits <- colnames(predresults[[1]])
tbyt <- list()

for (i in 1:length(gptraits)) {
  tmp <- lapply(predresults, function(x)
    data.frame(x[, i]))
  tmpout <- do.call(cbind, tmp)
  rownames(tmpout) <- rownames(predresults[[1]])
  colnames(tmpout) <- 1:ncol(tmpout)
  tbyt[[i]] <- tmpout
  names(tbyt)[i] <- gptraits[i]
}

nsubs <- ncol(tbyt[[1]])

p <- list()

# reshape the dataframes and plot them
for (i in 1:length(tbyt)){
  gptrait <- names(tbyt)[i]
  tmp <- tbyt[[i]]
  toplot <- melt(t(tmp))
  colnames(toplot) <- c("Rep", "Trait", "Accuracy")
  
  # get summary stats
  toplot %<>% 
    group_by(Trait) %>% 
    summarise(value = mean(Accuracy, na.rm = TRUE),
              se = sd(Accuracy) / sqrt(length(Accuracy)))
  
  p[[i]]<- ggplot(toplot) +
    geom_bar(aes(x=Trait, y=value), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar(aes(x=Trait, ymin=value-se, ymax=value+se), 
                  width=0.4, colour="black", alpha=0.9, size=0.5) +
    labs(title=gptrait)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
  
}


p[[5]]
