######################
# This scripts looks into the results of genomic selection 
# and produce all relevant plots
# Date: 01/05/2020
######################

# load needed packages
library(ggplot2)
library(ggpubr)
library(psych)
library(ggcorrplot)
library(patchwork)
library(corrplot)

#setting up a plotting function to be used thorughout
plotGP<-function(gpout, idx, predictors, predicted, cols=NA, gpname="", ymin=NULL, ymax=NULL){
  gpout<-gpout
  idx<-idx
  gp<-gpout[[idx]]
  
  gp<-gp[gp$trait %in% predicted, ]
  gp<-gp[gp$variable %in% predictors, ]
  
  #write input tables to file
  write.csv(gp, file=paste(names(gpout)[idx], "csv", sep="."), quote=F, row.names=F)
  
  if(is.na(gpname)==T){
    gpname<-names(gpout)[idx]
  }
  
  if(length(which(is.na(cols)))>0){
    cols<-gray.colors(length(predictors))
  }
  
  #get to plotting
  gp[,"value"]<-as.numeric(gp[,"value"])
  ggplot(gp, aes(x=trait, y=value, fill=variable)) + 
    geom_bar(stat="identity", position = "dodge", colour="black") + 
    scale_fill_manual(values=cols)+
    labs(fill = "variable") +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax))
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

######
#load all needed data 
#load BLUP data 
load("output/diversity.panel.BLUPs.Rdata")
mts<-blups.met
fms<-blups.farm
rownames(mts)<-sub("^ID_", "", mts[,1])
mts<-mts[,-1]
rownames(fms)<-sub("^ID_", "", fms[,1])
fms<-fms[,-1]

#load log abilities
load("output/log-abilities/all.traits.log-abilities.Rdata")

# load the results from the genomic selection
load("output/GP.output.Rdata")

#create lables with subscripts in a lazy way
gys<-expression('GY'[STATION])
gyf<-expression('GY'[FARM])
gysge<-expression('GY'[STATION]*'Ge')
gyshs<-expression('GY'[STATION]*'Hs')
gys12ge<-expression('GY'[STATION]*'2012.Ge')
gys12hs<-expression('GY'[STATION]*'2012.Hs')
gys13ge<-expression('GY'[STATION]*'2013.Ge')
gys13hs<-expression('GY'[STATION]*'2013.Hs')
gys12<-expression('GY'[STATION]*'2012')
gyf12<-expression('GY'[FARM]*'2012')
gys13<-expression('GY'[STATION]*'2013')
gyf13<-expression('GY'[FARM]*'2013')
gyf14<-expression('GY'[FARM]*'2014')
gyf15<-expression('GY'[FARM]*'2015')

oas<-expression('OA'[STATION])
oaf<-expression('OA'[FARM])
oasf<-expression('OA'[STATION]*"W")
oasm<-expression('OA'[STATION]*"M")
oasge<-expression('OA'[STATION]*'Ge')
oashs<-expression('OA'[STATION]*'Hs')
oas12ge<-expression('OA'[STATION]*'2012.Ge')
oas12hs<-expression('OA'[STATION]*'2012.Hs')
oas13ge<-expression('OA'[STATION]*'2013.Ge')
oas13hs<-expression('OA'[STATION]*'2013.Hs')
oas12<-expression('OA'[STATION]*'2012')
oaf12<-expression('OA'[FARM]*'2012')
oas13<-expression('OA'[STATION]*'2013')
oaf13<-expression('OA'[FARM]*'2013')
oaf14<-expression('OA'[FARM]*'2014')
oaf15<-expression('OA'[FARM]*'2015')

#fix genomic selection lables with a function
repnames<-function(x){
  for(i in 1:ncol(x)){
    tmp<-x[,i]
    #cent
    tmp<-sub("^GY","GY-Station",tmp)
    tmp<-sub("^OA","OA-Station",tmp)
    tmp<-sub("geregera","Ge",tmp)
    tmp<-sub("hagreselam","Hs",tmp)

    #decent
    tmp<-sub("d.OA", "OA-Farm", tmp)
    tmp<-sub("d.GY", "GY-Farm", tmp)
    tmp<-sub("^all.", "", tmp)
    x[,i]<-tmp
  }#for i
  x<-as.data.frame(x)
  return(x)
}#function

fixoutlist<-lapply(outlist,repnames)

####################
# produce plots

#correlation plot
#subset centralized data to samples in decentralized fields
mts<-mts[rownames(mts) %in% rownames(logab),]
fms<-fms[rownames(fms) %in% rownames(logab),]
stopifnot(all(rownames(mts)==rownames(fms)))
stopifnot(all(rownames(mts)==rownames(logab)))

#make a correlation df
cordf<-cbind(mts[,c("GY", "GY.hagreselam", "GY.geregera")], 
             fms[,c("OVERALL", "OVERALL.geregera", "OVERALL.hagreselam") ], 
             logab[,c("all.gy", "all.pow")])
colnames(cordf)<-c("GY-Station", "GY-Station.Hs", "GY-Station.Ge", "OA-Station", "OA-Station.Hs", "OA-Station.Ge", "GY-Farm", "OA-Farm")

cr<- corr.test(x=cordf, method ="spearman", 
          use="pairwise.complete")

corrplot(cr$r)

corplot<-ggcorrplot(cr$r, type=  "upper" , outline.col = "white", p.mat = cr$p, insig = "pch", show.diag=T)+
  geom_vline(xintercept=6.5, col="black", lty=2) + geom_hline(yintercept=6.5, col="black", lty=2)
   #geom_vline(xintercept=10.5, col="gray") + geom_hline(yintercept=10.5, col="gray")
corplot<- corplot + theme(legend.position="bottom") 
corplot<- corplot + scale_x_discrete(labels=c(gys, gyshs, gysge, oas, oashs, oasge, gyf, oaf))+ scale_y_discrete(labels=c(gys, gyshs, gysge, oas, oashs, oasge, gyf, oaf))
corplot

#ggsave("fig.s03.corplot.pdf", width = 9, height = 9)
#ggsave("fig.s03.corplot.png", width = 9, height = 9)

corplot<-corplot + labs(tag ="A")

#save tables
#write.csv(cr$r, file="Fig.2A.correlation.coefficients.csv", quote=F, row.names = F)
#write.csv(cr$p, file="Fig.2A.correlation.pvalues.csv", quote=F, row.names = F)

#station over station plot
idx<-1
#select traits to plot
predictors<-c("OA-Station",  "GY-Station.2012", "GY-Station.2012.Ge", "GY-Station.2012.Hs")
predicted<-c( "GY-Station.2013", "GY-Station.2013.Ge", "GY-Station.2013.Hs")
p1<-plotGP(gpout=fixoutlist, idx=idx, predictors=predictors, predicted=predicted)
p1<-addSmallLegend(p1) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p1<- p1 + labs(title = "Predicting stations", fill="Predictors", x ="Predicted", y = expression(tau)) #+ theme(aspect.ratio=0.5)
#fix lables
p1<- p1 + scale_x_discrete(labels=c(gys13, gys13ge, gys13hs)) + scale_fill_manual( labels = c(gys12, gys12ge, gys12hs, oas), values=gray.colors(length(predictors))) + theme(legend.text.align = 0)
p1

#ggsave("fig.s04.station.over.station.gp.pdf", width = 6, height = 4)
#ggsave("fig.s04.station.over.station.gp.png", width = 6, height = 4)

p1<-p1 + labs(tag="B")

#station over decentralized
idx<-2
predictors<-c("OA-Station",  "GY-Station")
predicted<-c("GY-Farm", "OA-Farm")
p2<-plotGP(gpout=fixoutlist, idx=idx, predictors=predictors, predicted=predicted, ymin=-0.2, ymax=0.3)
p2<-addSmallLegend(p2) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p2<- p2 + labs(title = "Predicting farms", fill="Predictors", x ="Predicted", y = expression(tau)) #+ theme(aspect.ratio=1)
#fix lables
p2<- p2 + scale_x_discrete(labels=c(gyf, oaf)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0)
p2 

#ggsave("fig.s06.station.over.farm.gp.pdf", width = 4, height = 4)
#ggsave("fig.s06.station.over.farm.gp.png", width = 4, height = 4)

p2<-p2 + labs(tag="C")


### produce a composite plot
comp<-(corplot | (p1 / (p2+ plot_spacer())))  
comp

#ggsave("fig.2.pdf", width = 12, height = 8)
#ggsave("fig.2.png", width = 12, height = 8)


###############SUPPLEMENTAL FIGURES########################

#station over station plot, by seasons
idx<-2
predictors<-c("OA-Station",  "GY-Station")
predicted<-c("2013.GY-Farm","2014.GY-Farm","2015.GY-Farm", "2013.OA-Farm","2014.OA-Farm","2015.OA-Farm")
p2season<-plotGP(gpout=fixoutlist, idx=idx, predictors=predictors, predicted=predicted)
p2season<-addSmallLegend(p2season) + geom_abline(intercept=0, slope=0, linetype="dashed", color="gray80")
p2season<- p2season + labs(title = "Predicting farms by season", fill="Predictors", x ="Predicted", y = expression(tau)) #+ theme(aspect.ratio=1)
#fix lables
p2season<- p2season + scale_x_discrete(labels=c(gyf13,gyf14,gyf15, oaf13,oaf14,oaf15)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0)
p2season 

#ggsave("station.over.farm.seasons.pdf", width = 6, height = 4)
#ggsave("station.over.farm.seasons.png", width = 6, height = 4)

#block validation GY
idx<-2
predictors<-c("OA-Station",  "GY-Station")
predicted<-c("2013.GY-Farm", "2014.GY-Farm", "2015.GY-Farm")
p3<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.2, ymax=0.3)
p3<- p3 + scale_x_discrete(labels=c(gyf13, gyf14, gyf15)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
  labs(title = "Block prediction of GY", fill="Predictors", x ="Predicted", y = expression(tau))

p3

#block validation OA
idx<-2
predictors<-c("OA-Station",  "GY-Station")
predicted<-c("2013.OA-Farm", "2014.OA-Farm", "2015.OA-Farm")
p4<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.2, ymax=0.3)
p4<-p4 + scale_x_discrete(labels=c(oaf13, oaf14, oaf15)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
    labs(title = "Block prediction of OA", fill="Predictors", x ="Predicted", y = expression(tau))

p4

blockgp <- p3|p4 & theme(legend.position = "right") 
blockgp + plot_layout(guides = "collect")

#ggsave("fig.s07.supplemental.block.validation.pdf", width = 6, height = 4)
#ggsave("fig.s07.supplemental.block.validation.png", width = 6, height = 4)

#man and women predictions
idx<-1
predictors<-c("OA-Station.W",  "OA-Station.M")
predicted<-c("GY-Station", "OA-Station")
p5<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.1, ymax=0.6)
p5<-p5 + scale_x_discrete(labels=c(gys, oas)) + scale_fill_manual(labels = c(oasm, oasf), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
      labs(title = "Station prediction by gender", fill="Predictors", x ="Predicted", y = expression(tau))

p5

idx<-2
predictors<-c("OA-Station.W",  "OA-Station.M")
predicted<-c("GY-Farm", "OA-Farm")
p6<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.1, ymax=0.6)
p6<-p6 + scale_x_discrete(labels=c(gyf, oaf)) + scale_fill_manual(labels = c(oasm, oasf), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
  labs(title = "Farm prediction by gender", fill="Predictors", x ="Predicted", y = expression(tau))
p6

mwgp <- p5|p6 & theme(legend.position = "right") 
mwgp + plot_layout(guides = "collect")

#ggsave("fig.s05.supplemental.men.women.pdf", width = 8, height = 4)
#ggsave("fig.s05.supplemental.men.women.png", width = 8, height = 4)

####hot and cold predictions
idx<-4
predictors<-c("GY-Station",  "OA-Station")
predicted<-c("2013.GY-Farm", "2014.GY-Farm", "2015.GY-Farm", "2013.OA-Farm", "2014.OA-Farm", "2015.OA-Farm")
p7<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.5, ymax=0.5)
p7<-p7 + scale_x_discrete(labels=c(gyf13, gyf14, gyf15, oaf13, oaf14, oaf15)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
  labs(title = "Cold adapted genotypes", fill="Predictors", x ="Predicted", y = expression(tau))
p7

idx<-5
predictors<-c("GY-Station",  "OA-Station")
predicted<-c("2013.GY-Farm", "2014.GY-Farm", "2015.GY-Farm", "2013.OA-Farm", "2014.OA-Farm", "2015.OA-Farm")
p8<-plotGP(gpout=fixoutlist,idx=idx, predictors=predictors, predicted=predicted, ymin=-0.5, ymax=0.5)
p8<-p8 + scale_x_discrete(labels=c(gyf13, gyf14, gyf15, oaf13, oaf14, oaf15)) + scale_fill_manual(labels = c(gys, oas), values=gray.colors(length(predictors)))+theme(legend.text.align = 0) +
  labs(title = "Warm adapted genotypes", fill="Predictors", x ="Predicted", y = expression(tau))
p8

tempgp <- p7|p8 & theme(legend.position = "right") 
tempgp + plot_layout(guides = "collect")

#ggsave("fig.s08.supplemental.temp.pdf", width = 8, height = 4)
#ggsave("fig.s08.supplemental.temp.png", width = 8, height = 4)


#########################################
#get numbers relative to predictions

names(fixoutlist)
####station on station
predictors<-c("OA-Station",  "GY-Station.2012", "GY-Station.2012.Ge", "GY-Station.2012.Hs")
predicted<-c( "GY-Station.2013", "GY-Station.2013.Ge", "GY-Station.2013.Hs")

tmp<-fixoutlist[["gp.cent"]]
tmp<-tmp[tmp[,1] %in% predicted,]
tmp<-tmp[tmp[,2] %in% predictors,]

range(tmp[,3])

#####means of decentralized breeding
tmp<-fixoutlist[["gp.decent"]]

#separate GY and OA as predictors
oa<-tmp[tmp[,2] == "OA-Station",]
oa<-oa[grep("^201[0-9].OA-Farm", oa[,1]),]
oa[,3]<-as.numeric(oa[,3])

gy<-tmp[tmp[,2] == "GY-Station",]
gy<-gy[grep("^201[0-9].GY-Farm", gy[,1]),]
gy[,3]<-as.numeric(gy[,3])

#GY as predictor
meanGYonGY<-mean(gy[,3])
sdGYonGY<-sd(gy[,3])

meanOAonOA<-mean(oa[,3])
sdOAonOA<-sd(oa[,3])

#OA as predictor
suboa<-tmp[grep("OA", tmp[,2]),]
meanOAonGY<-mean(as.numeric(suboa[grep("GY", suboa[,1]),3]))
sdOAonGY<-sd(as.numeric(suboa[grep("GY", suboa[,1]),3]))

meanOAonOA<-mean(as.numeric(suboa[grep("OA", suboa[,1]),3]))
sdOAonOA<-sd(as.numeric(suboa[grep("OA", suboa[,1]),3]))

#get means when using centralized OA as predictor
means<-aggregate(value ~ variable, FUN=mean, data=oa)
means
sds<-aggregate(value ~ variable, FUN=sd, data=oa)
sds

#get means when using centralized GY as predictor
means<-aggregate(value ~ variable, FUN=mean, data=gy)
means
sds<-aggregate(value ~ variable, FUN=sd, data=gy)
sds

#######################################
