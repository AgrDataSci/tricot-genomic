######################
#this script checks the environmental space covered by tested locations
##################

options(stringsAsFactors = F)
wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/4.env.space"
setwd(wd)

library(ggplot2)
library(plyr)
library(RCurl)
library(reshape2)
library(plotly)
library(plot3D)

#get station climate
url1<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices_station.csv")
stat<-read.csv(text=url1)
#get summary statistics
stat<-aggregate(. ~ location, data = stat, FUN = mean)
rownames(stat)<-stat[,1]
stat<-stat[,-1]
stat<-stat[,-c(58:ncol(stat))]

#get cs climate
url2<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices.csv")
farm<-read.csv(text=url2)
rownames(farm)<-farm$id
farm<-farm[,-c(49:ncol(farm))]

#merge the two dataframes only considering common columns
stat<-stat[,which(colnames(stat) %in% colnames(farm))]
farm<-farm[,which(colnames(farm) %in% colnames(stat))]
stopifnot(all(colnames(stat)==colnames(farm)))

#assemble the dataset and move to PCA
alldat<-data.frame(rbind(farm, stat))
colz<-rep("gray", nrow(alldat))
colz[(length(colz)-1):length(colz)]<-"red"

pca<-prcomp(alldat, scale = T)
pov<-pca$sdev^2/sum(pca$sdev^2)

plot(pca$x[,1:2], col=colz)

pctoplot<-data.frame(cbind(colz, pca$x[,1:3],alldat))
pctoplot$pch<- "circle"
pctoplot$pch[(length(pctoplot$pch)-1):length(pctoplot$pch)]<-"diamond-open"
pctoplot$cex<-"6"
pctoplot$cex[(length(pctoplot$cex)-1):length(pctoplot$cex)]<-"8"

#make a 3d pc plot with plotly
p <- plot_ly(pctoplot, x = ~PC1, y = ~PC2, z = ~PC3,
             marker = list(color = ~minDT_veg,
                           symbol = ~pch, 
                           size = 4, 
                           colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
             add_markers() %>%
              layout(scene = list(xaxis = list(title = paste('PC1', paste0(round(pov[1]*100,1), "%"))),
                      yaxis = list(title = paste('PC2', paste0(round(pov[2]*100,1), "%"))),
                      zaxis = list(title = paste('PC3', paste0(round(pov[3]*100,1), "%")))),
             annotations = list(
           x = 1.13,
           y = 1.05,
           text = 'minDT (°C)',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))
p


pctoplot2<-data.frame(cbind(colz, pca$x[,1:3],alldat))
pctoplot2$pch<- 20
pctoplot2$pch[(length(pctoplot2$pch)-1):length(pctoplot2$pch)]<-18
pctoplot2$cex<-1
pctoplot2$cex[(length(pctoplot2$cex)-1):length(pctoplot2$cex)]<-2

pdf("plot.3d.pdf")
plot3d<-scatter3D(pctoplot2[,"PC1"], pctoplot2[,"PC2"], pctoplot2[,"PC3"],
                  colvar=pctoplot2[,"minDT_veg"],  col = ramp.col(c("yellow","orange","red"),50),
                  pch = pctoplot2$pch, cex = pctoplot2$cex, 
                  bty = "g", colkey = list(side = 4, length = 0.5),
                  phi = 0, clab ="minDT_veg\n(°C)", 
                  bg="gray50") 
dev.off()

