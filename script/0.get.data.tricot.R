options(stringsAsFactors = F)

library(plyr)
library(RCurl)
library(lme4)
library(bestNormalize)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/0.get.data"
setwd(wd)

#get in raw data from GITHUB
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/durumwheat.csv")
csd<-read.csv(text=url)
csd<-csd[order(csd[,"genotype"]),]

#get a metric about getting first for each variety
bygen<-split(csd, csd[,"genotype"])
#set df as output
out.csd<-matrix(ncol=5, nrow=length(bygen))

for (i in 1:length(bygen)){
  tmp<-bygen[[i]]
  
  #get propotion of winning
  wins<-length(which(tmp[,"farmer_rank"]==1))
  odds<-wins/nrow(tmp)
  
  #get mean estimates
  gy<-mean(tmp[,"gy_gm"], na.rm=T)
  sl<-mean(tmp[,"mean_sl"], na.rm=T)
  nsps<-mean(tmp[,"mean_seed_no_spike"], na.rm=T)
  net<-mean(tmp[,"mean_tn"], na.rm=T)
  
  #bring everything in the output
  out.csd[i,]<-c(odds, gy, sl, nsps, net)
}

colnames(out.csd)<-c("win_prop", "mean_gy", "sl", "nsps", "net")
rownames(out.csd)<-unique(csd[,"genotype"])

out.csd<-as.data.frame(out.csd)

#get gy from a mixed model accounting for farmer and genotypes
mod<-lmer(gy_gm  ~  (1|genotype) + (1|farmer) + (1|year), data=csd)
idblup<-ranef(mod)
gy<-idblup$genotype
out.csd<-merge(out.csd, gy, by="row.names", all.x=T)
rownames(out.csd)<-out.csd[,1]
out.csd<-out.csd[,-1]
colnames(out.csd)[ncol(out.csd)]<-"gy_blup"

#now add in data from Kaue's models
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/log-abilities/log-abilities.csv")
rank<-read.csv(text=url)
rownames(rank)<-rank[,1]
rank<-rank[,-1]
out.csd<-merge(out.csd, rank, by="row.names", all.x=T)

row.names(out.csd)<-out.csd[,1]
out.csd<-out.csd[,-1]

#now make sure to normalize the data in the proper way
pvals<-apply(out.csd, 2, function(x) shapiro.test(x)$p)
pvals<-which(pvals<0.05)

#store it in a septarate column
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
save(out.csd, file="cs.data.estimates.Rdata")
