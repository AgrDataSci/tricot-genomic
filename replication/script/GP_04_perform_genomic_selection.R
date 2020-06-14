######################
# This scripts performs genomic selection using rrBLUP
# Date: 01/05/2020
######################

options(stringsAsFactors = F)

# load needed packages
library(ggplot2)
library(rrBLUP)
library(RCurl)
library(reshape2)
library(devtools)
#devtools::install_github("agrobioinfoservices/gosset", upgrade = "never")
library(gosset)
library(ggpubr)

####load log abilities 
load("output/log-abilities/all.traits.log-abilities.Rdata")
head(logab)
colnames(logab)<-sub(".pow", ".d.OA", colnames(logab))
colnames(logab)<-sub(".gy", ".d.GY", colnames(logab))
decent<-logab

####load SNP data
load("output/genoytpic.data.rrBLUP.Rdata")

####load BLUP data 
load("output/diversity.panel.BLUPs.Rdata")
mts<-blups.met
fms<-blups.farm

colnames(fms)<-sub("TILLER", "TL", colnames(fms))
colnames(fms)<-sub("OVERALL", "OA", colnames(fms))
colnames(fms)<-sub("SPIKE", "SP", colnames(fms))
colnames(fms)<-sub("EARLINESS", "ES", colnames(fms))
colnames(fms)<-sub(".F$", ".W", colnames(fms))
phenos<-merge(mts, fms, by="ID")
phenos[,1]<-sub("^ID_", "", phenos[,1])
rownames(phenos)<-phenos[,1]
phenos<-phenos[,-1]

cent<-phenos[,grep("OA|GY", colnames(phenos))]

#get correlations among and within datasets
cr<-cor(decent)
corrplot::corrplot(cr)
cr<-cor(cent)
corrplot::corrplot(cr)

################## GENOMIC SELECTION ####################

##############
#conduct GP of centralized over decentralized with sample overlap
#here we are using all 400 genotypes in centralized fields to predict performance of the 41 cs genotypes in decentralized fields


#test genotypes will be all
t.geno<-imputedcl
t.pheno<-cent
#validation genotypes will be those included in the decentralized
v.geno<-imputedcl[rownames(imputedcl) %in% rownames(decent),]
v.pheno<-decent

#run the loop by each phenotype
traitout<-list()

for (p in 1:ncol(t.pheno)){
  curp<-colnames(t.pheno)[p]
  print(curp)
  
  model<-mixed.solve(t.pheno[,p], Z=t.geno, K=NULL, SE= F, return.Hinv = F )
  effects<-model$u
  genoeffect <- v.geno %*% effects
  GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
  
  
  save(model, GEBV, file=paste0("./output/GEBV.", curp, ".station.on.farm.Rdata"))
  
  #use Pearson correlation to get accuracy
  #accuracy<-cor(GEBV, v.pheno, use="complete")
  
  #use Kendall Tau correlations from R/gosset to get accuracy
  #create a workaround to allow KT calcualtions
  prval<-rbind(GEBV, GEBV)
  toprval<-t(v.pheno)
  
  bytrait<-split(toprval, 1:nrow(toprval))
  #get a list for output
  bytrlist<-list()
  for(i in 1:length(bytrait)){
    tmp<-bytrait[[i]]
    tmptoprval<-rbind(tmp, tmp)
    bytrlist[[i]]<-kendallTau(prval, tmptoprval)
  }#for i
  
  accuracy<-data.frame(do.call(rbind, bytrlist))
  rownames(accuracy)<-rownames(toprval)
  
  traitout[[p]]<-accuracy[,1] #N effective can be caluclated here
  names(traitout[[p]])<- rownames(accuracy)
}  

gp.decent<-do.call(cbind, traitout)
colnames(gp.decent)<- colnames(t.pheno)

##############
#conduct GP of decentralized over centralized WITHOUT sample overlap
#here we are using the 41 cs genotypes in decentralized fields to predict performance of the non-overlapping 359 genotypes in centralized fields 

#test genotypes will be those NOT included in the decentralized
t.geno<-imputedcl[rownames(imputedcl) %in% rownames(decent),]
t.pheno<-decent
#validation genotypes will be those included in the decentralized
v.geno<-imputedcl[!rownames(imputedcl) %in% rownames(decent),]
v.pheno<-cent[!rownames(cent) %in% rownames(decent),]

#run the loop by each phenotype
traitout<-list()

for (p in 1:ncol(t.pheno)){
  curp<-colnames(t.pheno)[p]
  print(curp)
  
  model<-mixed.solve(t.pheno[,p], Z=t.geno, K=NULL, SE= F, return.Hinv = F )
  effects<-model$u
  genoeffect <- v.geno %*% effects
  GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
  
  save(model, GEBV, file=paste0("./output/GEBV.", curp, ".farm.on.station.Rdata"))
  
  #use Pearson correlation to get accuracy
  #accuracy<-cor(GEBV, v.pheno, use="complete")
  
  #use Kendall Tau correlations from R/gosset to get accuracy
  #create a workaround to allow KT calcualtions
  prval<-rbind(GEBV, GEBV)
  toprval<-t(v.pheno)
  
  bytrait<-split(toprval, 1:nrow(toprval))
  #get a list for output
  bytrlist<-list()
  for(i in 1:length(bytrait)){
    tmp<-bytrait[[i]]
    tmptoprval<-rbind(tmp, tmp)
    bytrlist[[i]]<-kendallTau(prval, tmptoprval)
  }#for i
  
  accuracy<-data.frame(do.call(rbind, bytrlist))
  rownames(accuracy)<-rownames(toprval)
  
  traitout[[p]]<-accuracy[,1]
  names(traitout[[p]])<- rownames(accuracy)
}  

gp.decent2cent<-do.call(cbind, traitout)
colnames(gp.decent2cent)<- colnames(t.pheno)

##############
#conduct GP of centralized over centralized on the samples that went in decentralized fields
#here we are using all 400 genotypes in centralized fields to predict performance of the 41 cs genotypes in centralized fields

#test genotypes will be all
t.geno<-imputedcl
t.pheno<-cent
#validation genotypes will be those included inm the CS, tested in the centralized fields
v.geno<-imputedcl[rownames(imputedcl) %in% rownames(decent),]
v.pheno<-cent[rownames(cent) %in% rownames(decent),]

#run the loop by each phenotype
traitout<-list()

for (p in 1:ncol(t.pheno)){
  curp<-colnames(t.pheno)[p]
  print(curp)
  
  model<-mixed.solve(t.pheno[,p], Z=t.geno, K=NULL, SE= F, return.Hinv = F )
  effects<-model$u
  genoeffect <- v.geno %*% effects
  GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
  
  save(model, GEBV, file=paste0("./output/GEBV.", curp, ".station.on.station.Rdata"))
 
  #use Pearson correlation to get accuracy
  #accuracy<-cor(GEBV, v.pheno, use="complete")
  
  #use Kendall Tau correlations from R/gosset to get accuracy
  #create a workaround to allow KT calcualtions
  prval<-rbind(GEBV, GEBV)
  toprval<-t(v.pheno)
  
  bytrait<-split(toprval, 1:nrow(toprval))
  #get a list for output
  bytrlist<-list()
  for(i in 1:length(bytrait)){
    tmp<-bytrait[[i]]
    tmptoprval<-rbind(tmp, tmp)
    bytrlist[[i]]<-kendallTau(prval, tmptoprval)
  }#for i
  
  accuracy<-data.frame(do.call(rbind, bytrlist))
  rownames(accuracy)<-rownames(toprval)
  
  traitout[[p]]<-accuracy[,1]
  names(traitout[[p]])<- rownames(accuracy)
}  

gp.cent<-do.call(cbind, traitout)
colnames(gp.cent)<- colnames(t.pheno)

##############
#conduct GP of centralized over farmer fields using only warm, adapted genotypes
#here we are using all 400 genotypes in centralized fields to predict performance of the warm apdapted cs genotypes in decentralized fields

load("output/phenotypes.by.adaptation.cs.samples.Rdata")

#sort out warm adapted samples
warms<-rownames(allph)[which(allph$adapt=="warm")]

#test genotypes will be those NOT included in the CS
t.geno<-imputedcl
t.pheno<-cent
#validation genotypes will be those included inm the CS, tested in the centralized fields
v.geno<-imputedcl[rownames(imputedcl) %in% warms,]
v.pheno<-decent[rownames(decent) %in% warms,]

#run the loop by each phenotype
traitout<-list()

for (p in 1:ncol(t.pheno)){
  curp<-colnames(t.pheno)[p]
  print(curp)
  
  model<-mixed.solve(t.pheno[,p], Z=t.geno, K=NULL, SE= F, return.Hinv = F )
  effects<-model$u
  genoeffect <- v.geno %*% effects
  GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
  
  #use Pearson correlation to get accuracy
  #accuracy<-cor(GEBV, v.pheno, use="complete")
  
  #use Kendall Tau correlations from R/gosset to get accuracy
  #create a workaround to allow KT calcualtions
  prval<-rbind(GEBV, GEBV)
  toprval<-t(v.pheno)
  
  bytrait<-split(toprval, 1:nrow(toprval))
  #get a list for output
  bytrlist<-list()
  for(i in 1:length(bytrait)){
    tmp<-bytrait[[i]]
    tmptoprval<-rbind(tmp, tmp)
    bytrlist[[i]]<-kendallTau(prval, tmptoprval)
  }#for i
  
  accuracy<-data.frame(do.call(rbind, bytrlist))
  rownames(accuracy)<-rownames(toprval)
  
  traitout[[p]]<-accuracy[,1]
  names(traitout[[p]])<- rownames(accuracy)
}  

gp.warm<-do.call(cbind, traitout)
colnames(gp.warm)<- colnames(t.pheno)


#####################
#Assemble the outputs and print to file

#for each of the gp outs, melt the outcome in the long format
outlist<-list(gp.cent, gp.decent, gp.decent2cent, gp.cold, gp.warm)
for (i in 1:length(outlist)){
  tmp<-data.frame(outlist[[i]])
  tmp$trait<-rownames(tmp)
  outlist[[i]]<-melt(tmp)
}#for i

names(outlist)<-c("gp.cent", "gp.decent", "gp.decent2cent", "gp.cold", "gp.warm")

#save everything as it is
save(outlist, file="output/GP.output.Rdata")


