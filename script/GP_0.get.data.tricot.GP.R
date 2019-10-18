# ..........................................
# ..........................................
# Combine genomic data with farmer metrics for genomic prediction
# Matteo Dell'Acqua
# ..........................................
# ..........................................

options(stringsAsFactors = F)

library(plyr)
library(RCurl)
library(lme4)
library(bestNormalize)
library(rrBLUP)

##################
#DIVERSITY PANEL DATA
##################

#load in the diversity panel datasets
load("data/diversity.panel.data.gp.Rdata")

#################
#SNP data

#set a function converting genotypic data in the proper format to run rrBLUP
convert<-function(x){
  tmp<-as.character(x)
  tmp[tmp=="N"]<-NA
  tmp[tmp  %in% c("R","Y","S","W","K","M") ]<-0
  uniqus<-sort(unique(tmp))
  uniqus<-uniqus[!uniqus %in% c("0", NA)]
  if(length(uniqus)>1){
    tmp[tmp==uniqus[1]]<- 1
    tmp[tmp==uniqus[2]]<- -1
  }    
  tmp[tmp==uniqus[1]]<-1
  #overwrite
  return(tmp)
}

#convert snps in -1,0,1 and NA
geno2<-apply(geno, 1, convert)
geno2<-data.frame(geno2)
geno2<-apply(geno2, 2, as.numeric)
rownames(geno2)<-colnames(geno)

#perform imputation with rrBLUP
imputed<-A.mat(geno2, max.missing=0.5, impute.method="mean", return.imputed=T)
failrate<-apply(imputed$imputed,2, function(x) length(which(is.na(x))))
todrop<-which(failrate>0)
if(length(todrop)>1){
  imputedcl<-imputed$imputed[,-which(failrate>0)]
} else {
  imputedcl<-imputed$imputed
}
rownames(imputedcl)<-colnames(geno)

output <- "processing/GP/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

save(imputedcl, file= paste0(output, "genoytpic.data.rrBLUP.Rdata"))

######################
#calculate BLUPs for metric traits

#get traits
met.traits<-colnames(met)[9:ncol(met)]

#calculate BLUP values for metric traits
metout<-list()
for (i in 1:length(met.traits)){
  trait<-met.traits[i]
  print(trait)
  formula<-as.formula(paste0(trait," ~ (1| ID) + (1 | LOCATION) + (1 | YEAR)  + (1|ID:YEAR:LOCATION)"))
  mod1<-lmer(formula, met)
  idblup<-ranef(mod1)
  metout[[i]]<-idblup$ID
}
metblup<-do.call(cbind, metout)
colnames(metblup)<-met.traits

#add BLUPs per location
metger<-met[met[,"LOCATION"]=="geregera",]
meths<-met[met[,"LOCATION"]=="hagreselam",]
listmet<-list(metger,meths)
names(listmet)<-c("geregera", "hagreselam")
listout<-list()
for (l in 1:length(listmet)){
  tmpmet<-listmet[[l]]
  metout<-list()
  for (i in 1:length(met.traits)){
    trait<-met.traits[i]
    print(trait)
    formula<-as.formula(paste0(trait," ~ (1| ID) + (1 | YEAR)"))
    mod1<-lmer(formula, tmpmet)
    idblup<-ranef(mod1)
    metout[[i]]<-idblup$ID
  }
  listout[[l]]<-do.call(cbind, metout)
  colnames(listout[[l]])<-paste(met.traits, names(listmet)[l], sep=".")
}

#assemble the dataset
for (i in 1:length(listout)){
  metblup<-merge(metblup, listout[[i]], by="row.names")
  rownames(metblup)<-metblup[,"Row.names"]
  metblup<-metblup[,-grep("Row.names", colnames(metblup))]
}

#add BLUPs per year
met12<-met[met[,"YEAR"]=="2012",]
met13<-met[met[,"YEAR"]=="2013",]
listmet<-list(met12,met13)
names(listmet)<-c("2012", "2013")
listout<-list()
for (l in 1:length(listmet)){
  tmpmet<-listmet[[l]]
  metout<-list()
  for (i in 1:length(met.traits)){
    trait<-met.traits[i]
    print(trait)
    formula<-as.formula(paste0(trait," ~ (1| ID) + (1 | LOCATION) "))
    mod1<-lmer(formula, tmpmet)
    idblup<-ranef(mod1)
    metout[[i]]<-idblup$ID
  }
  listout[[l]]<-do.call(cbind, metout)
  colnames(listout[[l]])<-paste(met.traits, names(listmet)[l], sep=".")
}

#assemble the dataset
for (i in 1:length(listout)){
  metblup<-merge(metblup, listout[[i]], by="row.names")
  rownames(metblup)<-metblup[,"Row.names"]
  metblup<-metblup[,-grep("Row.names", colnames(metblup))]
}

#add BLUPs for 2012 per location
ger12<-met[met[,"YEAR"]=="2012" & met[,"LOCATION"]=="geregera",]
hs12<-met[met[,"YEAR"]=="2012" & met[,"LOCATION"]=="hagreselam",]
listmet<-list(ger12,hs12)
names(listmet)<-c("geregera.2012", "hagreselam.2012")
listout<-list()
for (l in 1:length(listmet)){
  tmpmet<-listmet[[l]]
  metout<-list()
  for (i in 1:length(met.traits)){
    trait<-met.traits[i]
    print(trait)
    formula<-as.formula(paste0(trait," ~ (1|ID) + (1| REP) "))
    mod1<-lmer(formula, tmpmet)
    idblup<-ranef(mod1)
    metout[[i]]<-idblup$ID
  }
  listout[[l]]<-do.call(cbind, metout)
  colnames(listout[[l]])<-paste(met.traits, names(listmet)[l], sep=".")
}
#assemble the dataset
for (i in 1:length(listout)){
  metblup<-merge(metblup, listout[[i]], by="row.names")
  rownames(metblup)<-metblup[,"Row.names"]
  metblup<-metblup[,-grep("Row.names", colnames(metblup))]
}

#add BLUPs per 2013 per location
ger13<-met[met[,"YEAR"]=="2013" & met[,"LOCATION"]=="geregera",]
hs13<-met[met[,"YEAR"]=="2013" & met[,"LOCATION"]=="hagreselam",]
listmet<-list(ger13,hs13)
names(listmet)<-c("geregera.2013", "hagreselam.2013")
listout<-list()
for (l in 1:length(listmet)){
  tmpmet<-listmet[[l]]
  metout<-list()
  for (i in 1:length(met.traits)){
    trait<-met.traits[i]
    print(trait)
    formula<-as.formula(paste0(trait," ~ (1|ID) + (1| REP) "))
    mod1<-lmer(formula, tmpmet)
    idblup<-ranef(mod1)
    metout[[i]]<-idblup$ID
  }
  listout[[l]]<-do.call(cbind, metout)
  colnames(listout[[l]])<-paste(met.traits, names(listmet)[l], sep=".")
}

#assemble the dataset
for (i in 1:length(listout)){
  metblup<-merge(metblup, listout[[i]], by="row.names")
  rownames(metblup)<-metblup[,"Row.names"]
  metblup<-metblup[,-grep("Row.names", colnames(metblup))]
}

######################
#calculate BLUPs for farmer traits

#take out the infos before consolidating group scores
reducedfarm<-farm[,1:8]

#set number of groups per gender
ngroups<-3
#set traits
traits<-c("OV", "EAR", "SP", "TIL")
farm.traits<-c(paste0("overall", c(".m", ".w", "")),
               paste0("earliness", c(".m", ".w", "")),
               paste0("spike", c(".m", ".w", "")),
               paste0("tillering", c(".m", ".w", "")))

#get farmer means by trait, by gender
for (t in 1:length(traits)){
  print(traits[t])
  #men first
  m<-grep(paste0("^M[1-3]_", traits[t]), colnames(farm))
  #get means by farmer group
  listout<-list()
  for (i in 1:ngroups){
    tmpgend<-farm[,m]
    tmp<-tmpgend[,grep(paste0("M",i), colnames(tmpgend))]
    listout[[i]]<- rowMeans(tmp, na.rm=T)
  }
  men<-do.call(cbind, listout)
  colnames(men)<-paste0("M", 1:ngroups)
  #women now
  f<-grep(paste0("^F[1-3]_", traits[t]), colnames(farm))
  listout<-list()
  for (i in 1:ngroups){
    tmpgend<-farm[,f]
    tmp<-tmpgend[,grep(paste0("F",i), colnames(tmpgend))]
    listout[[i]]<- rowMeans(tmp, na.rm=T)
  }
  women<-do.call(cbind, listout)
  colnames(women)<-paste0("F", 1:ngroups)
  
  #assemble all in the final dataset
  reducedfarm<-cbind(reducedfarm, rowMeans(men, na.rm=T))
  reducedfarm<-cbind(reducedfarm, rowMeans(women, na.rm=T))
  reducedfarm<-cbind(reducedfarm, rowMeans(cbind(men, women), na.rm=T))
  
}#for t

#assign correct names
colnames(reducedfarm)[9:ncol(reducedfarm)]<-farm.traits

#now get BLUPs of farmer values                          
farmout<-list()
for (i in 1:length(farm.traits)){
  trait<-farm.traits[i]
  print(trait)
  formula<-as.formula(paste0(trait," ~ (1| ID) + (1 | LOCALITY) + (1 | REP)"))
  mod1<-lmer(formula, reducedfarm)
  idblup<-ranef(mod1)
  farmout[[i]]<-idblup$ID
}

#get out blup values
farmblup<-do.call(cbind, farmout)
colnames(farmblup)<-farm.traits

#add BLUPs per location
reducedfarmger<-reducedfarm[reducedfarm[,"LOCALITY"]=="geregera",]
reducedfarmhs<-reducedfarm[reducedfarm[,"LOCALITY"]=="hagreselam",]
listreducedfarm<-list(reducedfarmger,reducedfarmhs)
names(listreducedfarm)<-c("geregera", "hagreselam")
listout<-list()
for (l in 1:length(listreducedfarm)){
  tmpreducedfarm<-listreducedfarm[[l]]
  reducedfarmout<-list()
  for (i in 1:length(farm.traits)){
    trait<-farm.traits[i]
    print(trait)
    formula<-as.formula(paste0(trait," ~ (1| ID) + (1 | REP)"))
    mod1<-lmer(formula, tmpreducedfarm)
    idblup<-ranef(mod1)
    reducedfarmout[[i]]<-idblup$ID
  }
  
  tmpout<-reducedfarmout[[1]]
  for(j in 2:length(reducedfarmout)){
    tmpout<-merge(tmpout, reducedfarmout[[j]], by.x="row.names", by.y="row.names")
    rownames(tmpout)<-tmpout[,"Row.names"]
    tmpout<-tmpout[,-grep("Row.names", colnames(tmpout))]
  }
  colnames(tmpout)<-paste(farm.traits, names(listreducedfarm)[l], sep=".")
  listout[[l]]<-tmpout
}

#assemble the dataset
for (i in 1:length(listout)){
  farmblup<-merge(farmblup, listout[[i]], by="row.names")
  rownames(farmblup)<-farmblup[,"Row.names"]
  farmblup<-farmblup[,-grep("Row.names", colnames(farmblup))]
}

#bring data together and save it
phenos<-merge(metblup, farmblup, by="row.names")
rownames(phenos)<-phenos[,"Row.names"]
phenos<-phenos[,-grep("Row.names", colnames(phenos))]

save(phenos, file=paste0(output, "blup.diversity.panel.Rdata"))

##################
#CROWDSOURCING DATA
##################

#get in crowdsourcing data
#url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/durumwheat.csv")
#csd<-read.csv(text=url)
csd <- read.csv("data/durumwheat.csv")
csd<-csd[order(csd[,"genotype"]),]
csd <- csd[grepl("_D", csd$genotype), ]

#get infor about cs locations
#url1<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/data/environmental_indices.csv")
#csinfo<-read.csv(text=url1)

csinfo <- read.csv("data/environmental_indices.csv")

#categorize locations by climatic criteria emerging from the PL analysis
#namely, minDT_veg, splitting value 18.5?C
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
#url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/log-abilities/log-abilities.csv")
rank<-read.csv("output/log-abilities/log-abilities.csv")
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
save(out.csd, out.csd.regions, csids, file=paste0(output, "cs.data.estimates.Rdata"))

#in the end, produce an additive pairwise matrix from the SNP data
#subset genos to the cs materials
tokeep<-which(rownames(imputedcl) %in% csids)
genocs<-geno2[tokeep,]
#calculate the additive matrix and save it
genomatrix<-A.mat(genocs, max.missing=0.5, impute.method="mean", return.imputed=T)
additivemat<-genomatrix$A
save(additivemat, file= paste0(output, "additive.matrix.cs.rda"))
