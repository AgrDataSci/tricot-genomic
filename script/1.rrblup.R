######################
#this script performs GP with rrBLUP
#####################
options(stringsAsFactors = F)
library(ggplot2)
library(rrBLUP)
library(reshape2)
library(plyr)
library(RCurl)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/1.rrBLUP"
setwd(wd)

########  SNP DATA  ########
#get in genotypes and clean them
hmp<-read.delim("../../data/diversity.panel.svevo.positions.q10.hmp")
snp.pos<-hmp[,1:11]
geno<-hmp[,12:ncol(hmp)]
rownames(geno)<-hmp[,1]
colnames(geno)<-sub("^X", "", colnames(geno))

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
rownames(geno2)<-colnames(geno)
geno2<-apply(geno2, 2, as.numeric)

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

########  PHENOTYPIC DATA  ########
#get in BLUP values and clean them
phenos<-read.csv("../../data/Et_diversity_panel_ALL-TRAITS.BLUPS.csv")
phenos<-phenos[!is.na(phenos[,3]),]
rownames(phenos)<-phenos[,3]
phenos<-phenos[,6:ncol(phenos)]

#load in the crowdsourcing data produced by Kaue
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/exploring/probability_of_winning.csv")
csdata<-read.csv(text=url)
rownames(csdata)<-csdata[,1]
csdata<-csdata[,-1]

#make sure that the data is normal!!!
shapiro.test(csdata[,1])
shapiro.test(csdata[,2])

#add data produced locally
load("../0.get.data/rough.cs.data.estimates.Rdata")
#make it normally distributed
outrough<-apply(outrough,2,log)
joined.csdata<-merge(csdata, outrough, by="row.names")
row.names(joined.csdata)<-joined.csdata[,1]
joined.csdata<-joined.csdata[,-1]

#get crowsourcing IDs
csids<-rownames(joined.csdata)

#merge crowdsourcing with phenos
phenos<-merge(phenos, joined.csdata, by="row.names", all.x=T)
rownames(phenos)<-phenos[,1]
phenos<-phenos[,-1]

#remove all bread wheat samples!!!!!
phenos<-phenos[-grep("_B", rownames(phenos)), ]

########  DATA CLEANING  ########
#keep only samples that appear in all datasets
phenos<-phenos[which(rownames(phenos) %in% rownames(imputedcl)),]
imputedcl<-imputedcl[which(rownames(imputedcl) %in% rownames(phenos)),]

#sort everything right
phenos<-phenos[order(rownames(phenos)),]
imputedcl<-imputedcl[order(rownames(imputedcl)),]

#check that everything is OK
stopifnot(all(rownames(phenos)==rownames(imputedcl)))

#save to restart from here
save.image(file="GP.data.step1.Rdata")

######################################
#RESTART FROM HERE
options(stringsAsFactors = F)
library(rrBLUP)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/1.rrBLUP"
setwd(wd)

load("GP.data.step1.Rdata")

########  PERFORM GP  ########
#run GBLUP to perform genomic prediction
fraction<-.5
nruns<-100
tsize<-floor(length(csids)*fraction)

#generate samples
samps<-list()
for (i in 1:nruns){
  samps[[i]]<-sample(csids, tsize)
}
samps<-lapply(samps, sort)
if(length(which(duplicated(samps))>0)){
  samps[[which(duplicated(samps))]]<-NULL
}

#go on working on the genomic prediction
predresults<-list()

for (zz in 1:length(samps)){
  print(paste("iteration", zz, "of", length(samps)))
  testgeno<-samps[[zz]]
  traingeno<-setdiff(rownames(imputedcl), testgeno)

  #set training phenos and test phenos
  traintrait<-colnames(phenos)[grep("GY.201[23]$|SPIKE..2012$|OVERALL..2012$", colnames(phenos))]
  testtrait<-colnames(phenos)[97:ncol(phenos)]
  
  tpheno<-phenos[rownames(phenos) %in% traingeno,traintrait]
  tgeno<-imputedcl[rownames(imputedcl) %in% traingeno,]
  
  vpheno<-as.data.frame(phenos[rownames(phenos) %in% testgeno, testtrait])
  colnames(vpheno)<-testtrait
  rownames(vpheno)<-testgeno
  
  vgeno<-imputedcl[rownames(phenos) %in% testgeno, ]
  
  #sort traits to improve readability
  vpheno<-vpheno[,order(colnames(vpheno))]
  
  #run the loop by each phenotypes
  traitout<-list()
  for (p in 1:ncol(tpheno)){
    curp<-colnames(tpheno)[p]
    print(curp)
    
    trainpheno<-tpheno[,curp]
    
    model<-mixed.solve(trainpheno, Z=tgeno, K=NULL, SE= F, return.Hinv = F )
    effects<-model$u
    
    genoeffect <- vgeno %*% effects
    GEBV<-genoeffect[,1] + rep(model$beta, ncol(genoeffect))
    
    accuracy<-cor(GEBV, vpheno, use="complete")
    accuracy<-data.frame(t(accuracy))
    traitout[[p]]<-accuracy
    names(traitout[[p]])<-curp
  }
  traitoutdf<-do.call(cbind, traitout)
  colnames(traitoutdf)<-colnames(tpheno)
  predresults[[zz]]<-traitoutdf
  
}#for zz 

save(predresults, file="./GP.CSdata.Rdata")

######################################
#RESTART FROM HERE
options(stringsAsFactors = F)
library(ggplot2)
library(rrBLUP)
library(reshape2)
library(plyr)
library(RCurl)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/1.rrBLUP"
setwd(wd)

load("./GP.CSdata.Rdata")

#########produce some plotting########
gptraits<-colnames(predresults[[1]])
tbyt<-list()

for (i in 1:length(gptraits)){
  tmp<-lapply(predresults, function(x) data.frame(x[,i]))
  tmpout<-do.call(cbind, tmp)
  rownames(tmpout)<-rownames(predresults[[1]])
  colnames(tmpout)<-1:ncol(tmpout)
  tbyt[[i]]<-tmpout
  names(tbyt)[i]<-gptraits[i]
}

nsubs<-ncol(tbyt[[1]])

#reshape the dataframes and plot them
for (i in 1:length(tbyt)){
  gptrait<-names(tbyt)[i]
  tmp<-tbyt[[i]]
  toplot<-melt(t(tmp))
  colnames(toplot)<-c("Rep", "Trait", "Accuracy")
  
  #get summary stats
  toplot2<-ddply(toplot, ~ Trait, summarise, mean=mean(Accuracy, na.rm=T), sem=sd(Accuracy)/sqrt(length(Accuracy)))
  colnames(toplot2)<-c("Trait", "Accuracy", "SEM")
  
  png(paste0("./", gptrait, ".accuracy.", nsubs, ".repetitions.png"), width = 960, height = 480)
  barplot<-ggplot(toplot2) +
    geom_bar(aes(x=Trait, y=Accuracy), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar(aes(x=Trait, ymin=Accuracy-SEM, ymax=Accuracy+SEM), width=0.4, colour="black", alpha=0.9, size=0.5) +
    labs(title=gptrait)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(barplot)
  dev.off()
  
}#for t



############rubbish####################
outdf<-do.call(cbind, traitout)
colnames(outdf)<-colnames(tpheno)
outdf
meltdf<-melt(t(outdf), value.name = "Accuracy")
colnames(meltdf)<-c("Predictor", "Predicted", "Accuracy")

ggplot(meltdf, aes(x=Predicted,y=Accuracy,fill=factor(Predictor)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Predictor")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



