######################
#this script conducts analyses considering adaptation to DT_min_veg
#####################

options(stringsAsFactors = F)
wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/3.hot.cold"
setwd(wd)

library(ggplot2)
library(plyr)
library(RCurl)
library(reshape2)

#load diversity panel phentoypes
load("../0.get.data/blup.diversity.panel.Rdata")
#load crowdsourcing data
load("../0.get.data/cs.data.estimates.Rdata")
#load hot tolerance info
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/heat_tolerance/genotype_groups.csv")
hottol<-read.csv(text=url)
#set hot samples and cold samples
hot.samples<-hottol[which(hottol$group==1),1]
cold.samples<-hottol[which(hottol$group==2),1]

#add cold/hot info to cs data and dp data
allph<-merge(phenos, out.csd, by="row.names", all.y=T)
rownames(allph)<-allph[,1]
allph[,1]<-"warm"
allph[which(row.names(allph) %in% cold.samples),1]<-"cold"
colnames(allph)[1]<-"adapt"

allph$adapt<-as.factor(allph$adapt)

#sort out colds and hots
colds<-allph[which(allph$adapt=="cold"),2:ncol(allph)]
warms<-allph[which(allph$adapt=="warm"),2:ncol(allph)]


#check some correlations
cr<-cor(colds, use="complete.obs")
corrplot(cr)

#make some plots
boxplot(DF ~ adapt, data=allph)
t.test(DF ~ adapt, data=allph)

boxplot(GY ~ adapt, data=allph)
t.test(GY ~ adapt, data=allph)

boxplot(TGW ~ adapt, data=allph)
t.test(TGW ~ adapt, data=allph)

boxplot(overall ~ adapt, data=allph)
t.test(TGW ~ adapt, data=allph)

plot(allph$GY, allph$DF)

boxplot(GY_BLUP.cs ~ adapt, data=allph)
t.test(GY_BLUP.cs ~ adapt, data=allph)

#make violin plots
p <- ggplot(allph, aes(x=adapt, y=DB, fill=adapt)) + 
  geom_violin(trim=FALSE) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
  scale_fill_manual(values=c("#35E5FB", "#FA5B3C"))
p + theme(panel.background = element_rect(fill = 'white', colour = 'black'),axis.text.x = element_text(angle = 45, hjust = 1))
