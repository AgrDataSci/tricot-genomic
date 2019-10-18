# ..........................................
# ..........................................
# Genomic prediction with rrBLUP
# Matteo Dell'Acqua
# ..........................................
# ..........................................

options(stringsAsFactors = F)
library(ggplot2)
library(rrBLUP)
library(reshape2)
library(plyr)
library(RCurl)
library(corrplot)
library(psych)

#wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/1.rrBLUP"
#setwd(wd)

####load SNP data
load("processing/GP/genoytpic.data.rrBLUP.Rdata")

####load diversity panel phentoypes
#load("../0.get.data/blup.diversity.panel.Rdata")
mts<-read.delim("../../data/jesse_elaboration/diversity.panel.BLUPs.met.csv", sep=",")
fms<-read.delim("../../data/jesse_elaboration/diversity.panel.BLUPs.farm.csv", sep=",")
colnames(fms)<-sub("TILLER", "tillering", colnames(fms))
colnames(fms)<-sub("OVERALL", "overall", colnames(fms))
colnames(fms)<-sub("SPIKE", "spike", colnames(fms))
colnames(fms)<-sub("EARLINESS", "earliness", colnames(fms))
colnames(fms)<-sub(".F$", ".W$", colnames(fms))
phenos<-merge(mts, fms, by="ID")
phenos[,1]<-sub("^ID_", "", phenos[,1])
rownames(phenos)<-phenos[,1]
phenos<-phenos[,-1]

####load crowdsourcing data
load("../0.get.data/cs.data.estimates.Rdata")

####load hot tolerance info
url<-getURL("https://raw.githubusercontent.com/agrobioinfoservices/tricot-genomic/master/output/heat_tolerance/genotype_groups.csv")
hottol<-read.csv(text=url)

#start by bringing together crowdsourcing data and phenotypes
pheno<-merge(phenos, out.csd, by="row.names", all.x=T)
rownames(pheno)<-pheno[,1]
pheno<-pheno[,-1]

#make sure no bread wheat samples are included
todrop<-grep("_B", rownames(pheno))
if(length(todrop)>0){
  pheno<-pheno[-todrop, ]
}
dim(pheno)

#keep only samples that appear in all datasets
pheno<-pheno[which(rownames(pheno) %in% rownames(imputedcl)),]
imputedcl<-imputedcl[which(rownames(imputedcl) %in% rownames(pheno)),]

#sort everything right
pheno<-pheno[order(rownames(pheno)),]
imputedcl<-imputedcl[order(rownames(imputedcl)),]

#check that everything is OK
stopifnot(all(rownames(pheno)==rownames(imputedcl)))

######## MAKE SOME CORRELATIONS #########
csdonly<-pheno[rownames(pheno) %in% csids, ]
#select traits to correlate
forcor<-c("DF", "PH", "NET", "SPL", "SPS", "BM", "GY", "TGW", "overall", "spike",
          "win_prop.cs","SPL.cs", "SPS.cs", "NET.cs", "GY_BLUP.cs", "POW_rank", "POW_GY")

#cr<-cor(csdonly[,which(colnames(csdonly) %in% forcor)], method ="spearman", use="pairwise.complete")
cr<-corr.test(csdonly[,which(colnames(csdonly) %in% forcor)], method ="spearman", use="pairwise.complete")

#save correlation object and plot
save(cr, csdonly, forcor, file="correlation.data.Rdata")

pdf("correlations.pdf")
  corrplot(cr$r, p.mat = cr$p, insig = "n", "shade", tl.col="black")
dev.off()

########  PERFORM GP  ########
#run GBLUP to perform genomic prediction

#set trait space
dp.traits<-colnames(pheno)[1:126]
cs.traits<-colnames(pheno)[127:ncol(pheno)]
dp.2012<-colnames(pheno)[grep("over|spike|tiller|earl|2012", colnames(pheno))]
dp.2013<-colnames(pheno)[grep("2013", colnames(pheno))]
dp.farm.traits<-dp.traits[grep("over|spike|tiller|earl", dp.traits)]
dp.met.traits<-dp.traits[-grep("over|spike|tiller|earl", dp.traits)]

#put them in a list
tlist<-list(dp.traits, cs.traits, dp.2012, dp.2013, dp.farm.traits, dp.met.traits)
names(tlist)<-c("dp.traits", "cs.traits", "dp.2012", "dp.2013", "dp.farm.traits", "dp.met.traits")
  
#select tr1it space to be used for the prediction
idxtrain<-6
idxtest<-5

traintrait<-tlist[[idxtrain]]
testtrait<-tlist[[idxtest]]

#create name to be appended to your results
name<-paste(names(tlist)[idxtrain], "over", names(tlist)[idxtest], sep=".") 

#set the sample space from which to select test genotypes
all.samples<-rownames(imputedcl)
cs.samples<-csids
hot.samples<-hottol[which(hottol$group==1),1]
cold.samples<-hottol[which(hottol$group==2),1]

#put them in a list
slist<-list(all.samples, cs.samples, hot.samples, cold.samples)
names(slist)<-c("all.samples", "cs.samples", "hot.samples", "cold.samples")

#select trait space to be used as test for the prediction
idxsamp<-1
testgeno<-slist[[idxsamp]]

#set parameters of the GP routing
fraction<-0.7 #fraction of samples to be predicted (from the sample space)
tsize<-floor(length(testgeno)*fraction)
nruns<-100

#generate samples
samps<-list()
for (i in 1:nruns){
  samps[[i]]<-sample(testgeno, tsize)
}
#make sure samples are not repeated
samps<-lapply(samps, sort)
if(length(which(duplicated(samps))>0)){
  samps[[which(duplicated(samps))]]<-NULL
}

#go on working on the genomic prediction
predresults<-list()

#generate an output filename
outfilename<-paste0("./", name, ".", nruns, ".reps.", fraction, ".fraction.", names(slist)[[idxsamp]],".tested.jesse.data.GP.Rdata")
outfilename

#run the loop for number of iterations planned
for (zz in 1:length(samps)){
  print(paste("iteration", zz, "of", length(samps)))
  testgeno<-samps[[zz]]
  traingeno<-setdiff(rownames(imputedcl), testgeno)

  tpheno<-pheno[rownames(pheno) %in% traingeno,traintrait]
  tgeno<-imputedcl[rownames(imputedcl) %in% traingeno,]
  
  vpheno<-as.data.frame(pheno[rownames(pheno) %in% testgeno, testtrait])
  colnames(vpheno)<-testtrait
  rownames(vpheno)<-testgeno
  
  vgeno<-imputedcl[rownames(pheno) %in% testgeno, ]
  
  #sort traits to improve readability
  vpheno<-vpheno[,order(colnames(vpheno))]
  
  #run the loop by each phenotypes
  traitout<-list()
  for (p in 1:ncol(tpheno)){
    curp<-colnames(tpheno)[p]
    #print(curp)
    
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

save(predresults, file=outfilename)



######################################
#Make some plots

options(stringsAsFactors = F)
library(ggplot2)
library(rrBLUP)
library(reshape2)
library(plyr)
library(RCurl)

wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/1.rrBLUP"
setwd(wd)

#set a name to be appended to your results
#name<-"general"

load(paste0("./", name, ".", nruns, ".GP.CSdata.Rdata"))

#########produce some plotting########
gptraits<-colnames(predresults[[1]])
tbyt<-list()
listpred<-list()

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
listpreds<-list()

for (i in 1:length(tbyt)){
  gptrait<-names(tbyt)[i]
  tmp<-tbyt[[i]]
  toplot<-melt(t(tmp))
  colnames(toplot)<-c("Rep", "Trait", "Accuracy")
  
  #get summary stats
  toplot2<-ddply(toplot, ~ Trait, summarise, mean=mean(Accuracy, na.rm=T), sem=sd(Accuracy)/sqrt(length(Accuracy)))
  colnames(toplot2)<-c("Trait", "Accuracy", "SEM")
  toplot2<-cbind(rep(gptrait, nrow(toplot2)), toplot2)
  #assign to list
  listpreds[[i]]<-toplot2
  
  png(paste0("./",name, ".", gptrait, ".accuracy.", nsubs, ".repetitions.png"), width = 960, height = 480)
  barplot<-ggplot(toplot2) +
    geom_bar(aes(x=Trait, y=Accuracy), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar(aes(x=Trait, ymin=Accuracy-SEM, ymax=Accuracy+SEM), width=0.4, colour="black", alpha=0.9, size=0.5) +
    labs(title=gptrait)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(barplot)
  dev.off()
  
}#for t


#bring everything together and produce some plots
preddf<-do.call(rbind, listpreds)
colnames(preddf)<-c("Predictor", "Predicted", "Accuracy", "SEM")
rownames(preddf)<-NULL

#subset the dataframe prior plotting
predictors<-gptraits[grep("GY", gptraits)]
predicted<-c("mean_gy.norm", "pow_rank_bt", "nsps", "sl", "win_prop.norm")

#predicted<-c("pow_gy", "pow_rank_bt", "nsps", "sl", "win_prop.norm")

toplot<-preddf[preddf[,"Predicted"] %in% predicted, ]
toplot<-toplot[toplot[,"Predictor"] %in% predictors, ]

ggplot(toplot, aes(x=Predicted,y=Accuracy,fill=factor(Predictor)))+
  #scale_fill_manual(values=c(gray.colors(length(predictors))))+
  geom_bar(position=position_dodge(), aes(y=Accuracy), stat="identity") +
  geom_errorbar(aes(ymin=Accuracy-SEM, ymax=Accuracy+SEM), width=.2,position=position_dodge(.9)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),axis.text.x = element_text(angle = 45, hjust = 1)) 

