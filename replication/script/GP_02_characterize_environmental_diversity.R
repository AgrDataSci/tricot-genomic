######################
# This scripts derives environmental information in tested locations
# Date: 01/05/2020
######################

options(stringsAsFactors = F)

# load needed packages
library(ggplot2)
library(patchwork)
library(ggfortify)
library(RCurl)
library(rdist)
library(gtools)

#get station climate
url1<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices_station.csv")
stat<-read.csv(text=url1)
#get summary statistics
stat<-aggregate(. ~ location, data = stat, FUN = mean)
rownames(stat)<-stat[,1]
stat<-stat[,-1]
stat<-stat[,-c(58:ncol(stat))]
stat<-stat[,order(colnames(stat))]

#get farmer plots climate
url2<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices.csv")
farm<-read.csv(text=url2)
rownames(farm)<-farm$id
farm<-farm[,-c(49:ncol(farm))]
farm<-farm[,order(colnames(farm))]

#merge the two dataframes only considering common columns
stat<-stat[,which(colnames(stat) %in% colnames(farm))]
farm<-farm[,which(colnames(farm) %in% colnames(stat))]
stopifnot(all(colnames(stat)==colnames(farm)))

#assemble the dataset and move to PCA
alldat<-data.frame(rbind(farm, stat))
#subset to temperature variables
alldat<-alldat[,c("maxNT_rep","minDT_rep","minNT_sow2rep","minNT_veg", "DTR_gra","DTR_rep","DTR_sow2rep","DTR_veg")]

#########
# perform a multidimensional scaling and plot it

d <- dist(alldat) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

mds1<-data.frame(loc=rownames(alldat), MDS_1=fit$points[,1], MDS_2=fit$points[,2])
mds1$Location<-rep("Farmer field", nrow(mds1))
mds1$Location[1166:1167]<-c("Geregera", "Hagreselam")

m1<-ggplot(mds1, aes(x=MDS_1))+geom_histogram(fill="gray60")+
       geom_point(aes(color = Location, x=MDS_1, y=0),size = 4, data = mds1[mds1$Location != "Farmer field", ]) +
       theme_bw()

m2<-ggplot(mds1, aes(x=MDS_2))+geom_histogram(fill="gray60")+
  geom_point(aes(color = Location, x=MDS_2, y=0),size = 4, data = mds1[mds1$Location != "Farmer field", ]) +
  theme_bw()

mdsplot <- m1|m2 & theme(legend.position = "right") 
mdsplot + plot_layout(guides = "collect")

#ggsave("fig.s09.supplemental.mds.histogram.pdf", width = 9, height = 6)
#ggsave("fig.s09.supplemental.mds.histogram.png", width = 9, height = 6)


#########
# perform a PCA and plot it

pca<-prcomp(alldat, scale = T)
pov<-pca$sdev^2/sum(pca$sdev^2)

dftoplot<-data.frame(Location=mds1$Location, minDT_rep=alldat$minDT_rep, PC1=pca$x[,1], PC2=pca$x[,2])

pcplt<-autoplot(pca, scale = 0, data=dftoplot, colour="maxNT_rep") + 
        scale_color_gradient2(midpoint=12.5, low = "yellow", mid = "orange", high = "red", space = "Lab" )+
        geom_point(data=dftoplot[1166:1167,], aes(x=PC1, y=PC2, alpha=""), size=3,  fill = "gray", pch = 21)+
        labs(title="Variation in temperature", colour="max NT (Â°C)") +
         scale_alpha_manual(values = 1) +
        labs(alpha = "Stations")+
        theme_bw()  
pcplt

#ggsave("fig.s10.supplemental.pca.pdf", width = 6, height = 4)
#ggsave("fig.s10.supplemental.pca.png", width = 6, height = 4)

###########
# get warm adapted/cold adapted accessions

#get in decentralized data
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/durumwheat.csv")
csd<-read.csv(text=url)
csd<-csd[order(csd[,"genotype"]),]

#get infor about cs locations
url1<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices.csv")
csinfo<-read.csv(text=url1)

#categorize locations by climatic criteria emerging from the PL analysis
#namely, minDT_veg, splitting value 18.5°C
splitting<-18.5
hotsites<-csinfo$id[which(csinfo$minDT_veg>=splitting)]
coldsites<-csinfo$id[which(csinfo$minDT_veg<splitting)]

#divide the cs dataset in hot and cold
hotcsd<-csd[csd$id %in% hotsites,]
coldcsd<-csd[csd$id %in% coldsites,]

#get variable summary for each
csdlist<-list(csd, hotcsd, coldcsd)
names(csdlist)<-c("cs", "cs.hot", "cs.cold")
csdout<-list()

for(dat in 1:length(csdlist)){
  tmpcsd<-csdlist[[dat]]
  
  #get a metric about getting first for each variety
  bygen<-split(tmpcsd, tmpcsd[,"genotype"])
  
  #set df as output
  out.csd<-matrix(ncol=6, nrow=length(bygen))
  
  for (i in 1:length(bygen)){
    tmp<-bygen[[i]]
    tmpreg<-unique(tmp$region)
    if(length(tmpreg)>1){
      tmpreg<-paste(tmpreg, collapse="_")
    }
    
    #get propotion of winning
    wins<-length(which(tmp[,"farmer_rank"]==1))
    odds<-wins/nrow(tmp)
    
    #get mean estimates
    gy<-mean(tmp[,"gy_gm"], na.rm=T)
    sl<-mean(tmp[,"mean_sl"], na.rm=T)
    nsps<-mean(tmp[,"mean_seed_no_spike"], na.rm=T)
    net<-mean(tmp[,"mean_tn"], na.rm=T)
    
    #bring everything in the output
    out.csd[i,]<-c(tmpreg, odds, gy, sl, nsps, net)
  }
  
  colnames(out.csd)<-paste(c("region","win_prop", "GY.mean", "SPL", "SPS", "NET"), names(csdlist[dat]), sep=".")
  rownames(out.csd)<-unique(csd[,"genotype"])
  
  out.csd<-as.data.frame(out.csd)
  out.csd[,2:ncol(out.csd)]<-apply(out.csd[,2:ncol(out.csd)], 2, as.numeric)
  
  #get gy from a mixed model accounting for farmer and genotypes
  mod<-lmer(gy_gm  ~  (1|genotype) + (1|farmer) + (1|year), data=tmpcsd)
  idblup<-ranef(mod)
  gy<-idblup$genotype
  out.csd<-merge(out.csd, gy, by="row.names", all.x=T)
  rownames(out.csd)<-out.csd[,1]
  out.csd<-out.csd[,-1]
  colnames(out.csd)[ncol(out.csd)]<-paste("GY_BLUP", names(csdlist[dat]), sep=".")
  
  csdout[[dat]]<-out.csd
} #for dat    

#assemble in one object
out.csd<-do.call(cbind, csdout)
out.csd<-out.csd[,-grep("region.cs.", colnames(out.csd))]

#now add in data from PL models
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/log-abilities/log-abilities.csv")
rank<-read.csv(text=url)
rownames(rank)<-rank[,1]
rank<-rank[,-1]
colnames(rank)<-c("POW_rank", "POW_rank.bt", "POW_GY", "POW_GY.bt")

out.csd<-merge(out.csd, rank, by="row.names", all.x=T)

row.names(out.csd)<-out.csd[,1]
out.csd.regions<-data.frame(out.csd[,1:2])
out.csd<-out.csd[,-c(1,2)]

#get crowsourcing IDs
csids<-rownames(out.csd)
#make sure everything is numeric and bring back rownames
out.csd<-apply(out.csd, 2, as.numeric)
rownames(out.csd)<-csids

#now make sure to normalize the data in the proper way
pvals<-apply(out.csd, 2, function(x) shapiro.test(x)$p)
pvals<-which(pvals<0.05)

#store it in a separate column
out.csdn<-out.csd[,pvals]
idx<-1
for (i in pvals){
  x<- bestNormalize(out.csd[,i], allow_lambert_s = TRUE)
  out.csdn[,idx] <- predict(x)
  idx<-idx+1
}

colnames(out.csdn)<-sub("$", ".norm", colnames(out.csdn))

#put everything together and save it
out.csd<-cbind(out.csd, out.csdn)
save(out.csd, out.csd.regions, csids, file="output/cs.data.estimates.Rdata")

#load diversity panel phentoypes
load("output/diversity.panel.BLUPs.Rdata")
#load crowdsourcing data
load("output/cs.data.estimates.Rdata")
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

save(allph, file="output/phenotypes.by.adaptation.cs.samples.Rdata")

###########
# derive quantiles for environmental distances

#get  climate of centralized stations
url1<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices_station.csv")
cent.clim<-read.csv(text=url1)
#get summary of centralized climate
cent.clim<-aggregate(. ~ location, data = cent.clim, FUN = mean)
rownames(cent.clim)<-cent.clim[,1]
cent.clim<-cent.clim[,-1]

#get climate of decentralized trials
url2<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices.csv")
decent.clim<-read.csv(text=url2)
rownames(decent.clim)<-decent.clim$id

#merge the two dataframes only considering common columns
cent.clim<-cent.clim[,which(colnames(cent.clim) %in% colnames(decent.clim))]
decent.clim<-decent.clim[,which(colnames(decent.clim) %in% colnames(cent.clim))]
cent.clim<-cent.clim[,order(colnames(cent.clim))]
decent.clim<-decent.clim[,order(colnames(decent.clim))]

stopifnot(all(colnames(cent.clim)==colnames(decent.clim)))

#assemble the dataset and compute multivariate statistics
alldat<-data.frame(rbind(decent.clim, cent.clim))

#get pairwise distances using selected environemtnal info
alldat<-alldat[,c("maxNT_rep","minDT_rep","minNT_sow2rep","minNT_veg", "DTR_gra","DTR_rep","DTR_sow2rep","DTR_veg")]

d<-pdist(alldat, "manhattan")
colnames(d)<-rownames(alldat)
rownames(d)<-colnames(d)
dvsc<-d[,c("geregera", "hagreselam")]

#get an index of environmental distance
gerqt<-as.numeric(quantcut(dvsc[,"geregera"], q=4)) #using quartiles of the distribution
names(gerqt)<-rownames(dvsc)

hsqt<-as.numeric(quantcut(dvsc[,"hagreselam"], q=4))
names(hsqt)<-rownames(dvsc)

envqt<-data.frame(id=rownames(dvsc), ger.qt=gerqt, hs.qt=hsqt)
save(envqt, file="output/environmental.quantiles.by.centralized.loc.Rdata")
