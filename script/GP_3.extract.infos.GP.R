######################
#this script looks into predictions and extracts informations
#####################

options(stringsAsFactors = F)
wd<-"C:/Users/admin/OneDrive - Scuola Superiore Sant'Anna/projects/tricot.genomics/analyses/2.extract.info"
setwd(wd)

library(ggplot2)
library(plyr)
library(reshape2)
library(ggpubr)
library(psych)
library(ggcorrplot)

#start getting in heritability values
h2<-read.delim("../../data/jesse_elaboration/diversity.panel.h2.met.csv", sep=",")
h2$code<-paste(h2[,"trait"], h2[,"year"], h2[,"location"], sep="_")

#list output files in rrBLUP folder
filen<-list.files("../1.rrBLUP/")
filen<-filen[grep("GP.Rdata", filen)]

#get only those with Jesse's BLUPs
filen<-filen[grep("jesse", filen)]

#load the results of GP according to listed files
rslist<-list()
for(i in 1:length(filen)){
  load(paste0("../1.rrBLUP/",filen[i]))
  rslist[[i]]<-predresults
}#for i
names(rslist)<-filen
names(rslist)

#make an output list
outlist<-list()

#create output vector holding ALL trait nanmes
alltraitnames<-c()

#summarize results by trait by prediction
for(pl in 1:length(rslist)){
  
  tmppred<-rslist[[pl]]
  
  gptraits<-colnames(tmppred[[1]])
  
  alltraitnames<-c(alltraitnames, gptraits)
  
  tbyt<-list()
  
  for (i in 1:length(gptraits)){
    tmp<-lapply(tmppred, function(x) data.frame(x[,i]))
    tmpout<-do.call(cbind, tmp)
    rownames(tmpout)<-rownames(tmppred[[1]])
    colnames(tmpout)<-1:ncol(tmpout)
    tbyt[[i]]<-tmpout
    names(tbyt)[i]<-gptraits[i]
  }
  
  outlist[[pl]]<-tbyt
  
}#for pl
names(outlist)<-names(rslist)

#get unique trait names
alltraitnames<-sort(unique(alltraitnames))

#set up colors matching trait names
coldf<-data.frame(alltraitnames, col=NA)
#set color scales
grs  <- colorRampPalette(c("green1", "green4"))
grs2 <- colorRampPalette(c("olivedrab1", "olivedrab4"))
grs3  <- colorRampPalette(c("seagreen1", "seagreen4"))

bls <- colorRampPalette(c("blue1", "blue4"))
bls2 <- colorRampPalette(c("slateblue1", "slateblue4"))
bls3 <- colorRampPalette(c("lightskyblue1", "lightskyblue4"))
bls4 <- colorRampPalette(c("royalblue1", "royalblue4"))

rds <- colorRampPalette(c("red1", "red4"))
rds2<-  colorRampPalette(c("orangered1", "orangered4"))
rds3<- colorRampPalette(c("firebrick1", "firebrick4"))
rds4 <- colorRampPalette(c("tomato1", "tomato4"))

prp<- colorRampPalette(c("purple1", "purple4"))
prp2<- colorRampPalette(c("mediumorchid1", "mediumorchid4"))
prp3<- colorRampPalette(c("violetred1", "violetred4"))
prp4<- colorRampPalette(c("deeppink1", "deeppink4"))

yw<-colorRampPalette(c("yellow1", "yellow4"))
yw2<-colorRampPalette(c("goldenrod1", "goldenrod4"))
yw3<-colorRampPalette(c("lightsalmon1", "lightsalmon2"))


#order in file is
#BM, DB, DF, DM, earl, 
#GY, NET, over, PH, spike, 
#SPL, SPS, TGW, till

allcol<-c(rds(9), grs2(9), grs3(9), grs(9), bls(9), 
          rds2(9), yw(9), bls2(9), prp(9), bls3(9), 
          prp3(9), prp4(9), rds4(9), bls4(9))

  
coldf[,2]<-allcol
write.table(coldf, file="coldf.txt", sep="\t", quote=F, row.names = F)

#setting up a plotting function to be used thorughout
plotGP<-function(gpout, predictors, predicted, cols=NA, gpname="", ymin=-0.5, ymax=0.5){
  gp<-outlist[[gpout]]
  
  if(is.na(gpname)==T){
    gpname<-names(outlist)[gpout]
  }
  
  trin<-names(gp)
  trout<-rownames(gp[[1]])
  
  if(length(which(is.na(cols)))>0){
    cols<-gray.colors(length(predictors))
  }
  
  #get to plotting
  toplot<-list()
  for (i in 1:length(predictors)){
    tmppred<-trin[which(trin==predictors[i])]
    predictdf<-gp[[tmppred]]
    tmpdf<-predictdf[predicted,]
    tmpdf<-melt(t(tmpdf))
    colnames(tmpdf)<-c("Rep", "Trait", "Accuracy")
    
    #get summary stats
    tmpdf2<-ddply(tmpdf, ~ Trait, summarise, mean=mean(Accuracy, na.rm=T), sem=sd(Accuracy)/sqrt(length(Accuracy)))
    colnames(tmpdf2)<-c("Trait", "Accuracy", "SEM")
    tmpdf2$Predictor<-predictors[i]
    
    toplot[[i]]<-tmpdf2
  }#for i
  names(toplot)<-predictors
  
  #join dataframes
  if(length(toplot)>1){
    df<-do.call(rbind, toplot)
  } else {
    df<-toplot[[i]]
  }
  
  #write down the dataframe
  save(df, file=paste0("barplot.df.",names(outlist)[gpout]))
  
  #make the plot fixing the dodging
  pos <- position_dodge(.9)
  lim <- aes(ymin=Accuracy-SEM, ymax=Accuracy+SEM)
  
  barplot<-ggplot(df,aes(x=Trait, y=Accuracy, fill=factor(Predictor, levels=predictors))) + 
    scale_fill_manual(values=cols)+
    ggtitle(gpname) +
    labs(fill = "Predictors") +
    geom_bar(stat="identity", position = "dodge", colour="black") + 
    geom_bar(stat="identity",position=pos) + 
    geom_errorbar(lim, position=pos, width=.5) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax))
  print(barplot)
}### end of function

#include function to reduce legend size
addSmallLegend <- function(myPlot, pointSize = 0.9, textSize = 9, spaceLegend = 0.8) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

#######################
#make plots
######

#plot correlations
load("../1.rrBLUP/correlation.data.Rdata")
cr<-corr.test(csdonly[,which(colnames(csdonly) %in% forcor)], method ="spearman", 
              use="pairwise.complete")

corplot<-ggcorrplot(cr$r, outline.col = "white", p.mat = cr$p, insig = "blank")+
       #geom_vline(xintercept=c(0.5, 17.5)) + geom_hline(yintercept=c(0.5, 17.5))+
         geom_vline(xintercept=10.5, col="gray") + geom_hline(yintercept=10.5, col="gray")
corplot

#plot h2s
#make a plot siding h2 in GY and in farmer scores
h2sub<-h2[grep("GY",h2[,5]),]
h2sub$Gender<-"-"
h2sub$Type<-"MT"

#artifically incorporate farmers' h2 while they are calcualted
h2f<-data.frame(trait=rep("Overall", 9), year=rep("2012", 9),
                          location= c(rep("ALL", 3), rep("Geregera", 3), rep("Hagreselam", 3)),
                          h2=rep(NA,9))
h2f$code<-paste(h2f[,"trait"], h2f[,"year"], h2f[,"location"], sep="_")
h2f$Gender=rep(c("Both", "Men", "Women"),3)
h2f$Type<-rep("PVS", 9)
  
h2f[,4]<-c(0.52,0.36,0.625,0.635,0.47,0.56,0.32,0.18,0.525)
h2plot<-rbind(h2sub,h2f)
h2plot<-data.frame(apply(h2plot,2,function(x) sub("ALL", "Combined", x)))
h2plot<-data.frame(apply(h2plot,2,function(x) sub("geregera", "Geregera", x)))
h2plot<-data.frame(apply(h2plot,2,function(x) sub("hagreselam", "Hagreselam", x)))

h2plot[,4]<-as.numeric(h2plot[,4])

hplot<-ggplot(data=h2plot, aes(x=code, y=h2, fill=Gender)) +
       geom_bar(stat="identity", position="dodge", color="black")+  
       facet_wrap( ~location, ncol=1.,scales = "free") +
       theme(panel.background = element_rect(fill = 'white', colour = 'black'), 
            legend.position ="bottom", axis.title.y = element_blank(), legend.box="horizontal")+
       coord_flip() + 
        scale_fill_manual(values=c("gray90",'darkgrey','skyblue','pink'))
hplot<-addSmallLegend(hplot)
hplot

#plot GP results
names(outlist)

#########PLOT 1
#station over station plot
idx<-1
#select traits to plot
predictors<-c("GY.2012", "DF.2012", "TGW.2012","overall", "spike", "earliness")
predicted<-c("GY.2013" , "DF.2013","TGW.2013" , "BM.2013", "SPS.2013", "SPL.2013")
#set color scale
colz<- mapvalues(predictors, from=coldf[,1], to=coldf[,2], warn_missing=F)
p1<-plotGP(gpout=idx, predictors=predictors, predicted=predicted, cols=colz)
p1<-addSmallLegend(p1)

#########PLOT 2
#station over CS plot
idx<-6
#select traits to plot
predictors<-c("GY", "GY.hagreselam", "GY.geregera", "overall", "spike")
predicted<-c("GY_BLUP.cs", "win_prop.cs", "POW_GY", "POW_rank")
#set color scale
colz<- mapvalues(predictors, from=coldf[,1], to=coldf[,2], warn_missing=F)
p2<-plotGP(gpout=idx, predictors=predictors, predicted=predicted,  cols=colz) 
p2<-addSmallLegend(p2)

#########PLOT 3
#plot prediction ability in COLD varieties
idx<-8
#select traits to plot
predictors<-c("GY", "overall", "spike")
predicted<-c("GY_BLUP.cs", "GY_BLUP.cs.hot", "GY_BLUP.cs.cold", "win_prop.cs", "win_prop.cs.hot", "win_prop.cs.cold")
#set color scale
colz<- mapvalues(predictors, from=coldf[,1], to=coldf[,2], warn_missing=F)
p3<-plotGP(gpout=idx, predictors=predictors, predicted=predicted,  cols=colz, ymin=-0.1, ymax=0.8)
p3<-addSmallLegend(p3)

#########PLOT 4
#plot prediction ability in HOT varieties
idx<-10
#select traits to plot
predictors<-c("GY", "overall", "spike")
predicted<-c("GY_BLUP.cs", "GY_BLUP.cs.hot", "GY_BLUP.cs.cold", "win_prop.cs", "win_prop.cs.hot", "win_prop.cs.cold")
#set color scale
colz<- mapvalues(predictors, from=coldf[,1], to=coldf[,2], warn_missing=F)
p4<-plotGP(gpout=idx, predictors=predictors, predicted=predicted,  cols=colz, ymin=-0.5, ymax=0.5)
p4<-addSmallLegend(p4)

######Assemble all plots
#make a general plot
corplot.space<-corplot + theme(plot.margin = unit(c(1,1,1,3), "lines"))
upplot<-ggarrange(hplot, corplot.space, ncol = 2, heights=c(1, 0.8), widths =c(0.5, 0.8),
                  labels = c("A", "B"))
upplot

bottomplot<-ggarrange(p1 + rremove("xlab") + ggtitle("Station across seasons"), 
                  p2 + rremove("xylab") + ggtitle("Station over CS"), 
                  p3 + rremove("legend.title") + ggtitle("Below-threshold DTmin"), 
                  p4 + rremove("ylab") + ggtitle("Above-threshold DTmin"),  
                  labels = c("C", "D", "E", "F"), common.legend = F, legend="right")
bottomplot

allplot<-ggarrange(upplot,
                   bottomplot,
                   nrow=2)

png("composite.plot.v2.png", width=800, height = 1000)
  print(allplot)
dev.off()



#pdf("composite.plot.pdf", height=8, width=16)
#  print(allplot)
#dev.off()

