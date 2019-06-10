options(stringsAsFactors = F)

library(plyr)
library(RCurl)
library(lme4)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/0.get.data"
setwd(wd)

#get in row data from GITHUB
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/durumwheat.csv")
csd<-read.csv(text=url)
csd<-csd[order(csd[,"genotype"]),]

#try to drop samples from OROMIA
csd<-csd[which(csd[,"region"] %in% c("Amhara", "Tigray")),]

#get a metric about getting first for each variety
bygen<-split(csd, csd[,"genotype"])
#set df as output
outrough<-matrix(ncol=2, nrow=length(bygen))

for (i in 1:length(bygen)){
  tmp<-bygen[[i]]
  
  #get propotion of winning
  wins<-length(which(tmp[,"farmer_rank"]==1))
  odds<-wins/nrow(tmp)
  
  #get mean estimate
  gy<-mean(tmp[,"gy_gm"], na.rm=T)

  #bring everything in the output
  outrough[i,]<-c(odds, gy)
}

colnames(outrough)<-c("win_prop", "mean_gy")
rownames(outrough)<-unique(csd[,"genotype"])

outrough<-as.data.frame(outrough)

#get gy from a mixed model accounting for farmer and genotypes
mod<-lmer(gy_gm  ~  (1|genotype) + (1|farmer_no) + (1|year), data=csd)
idblup<-ranef(mod)
gy<-idblup$genotype

outrough<-merge(outrough, gy, by="row.names", all.x=T)
rownames(outrough)<-outrough[,1]
outrough<-outrough[,-1]
colnames(outrough)[3]<-"gy_blup"

save(outrough, file="rough.cs.data.estimates.Rdata")
