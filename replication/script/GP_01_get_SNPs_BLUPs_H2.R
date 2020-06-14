######################
# This scripts derives best linear unbiased predictions (BLUPs)
# and heritability estimates for the data collected in centralized stations
# Date: 01/05/2020
######################

options(max.print=10000)

# load needed packages
library(sm)
library(vioplot)
library(corrplot)
library(tidyr)
library(asreml)
library(bestNormalize)
library(rrBLUP)

#crate output folder
output <- "output"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# state functions to derive BLUPs, BLUEs, H2
etBLUPs = function(m){
    res = data.frame(coef(m)$rand)
    res = res[grepl("ID",row.names(res)) & !grepl(":",row.names(res)) ,  , drop=FALSE]
    ##res = res[grepl("GY:ID_",row.names(res)), , drop=FALSE]
    return(res)
}

getBLUEs = function(m){
    res = data.frame(coef(m)$fix)
    res = res[grepl("ID",row.names(res)) & !grepl(":",row.names(res)) ,  , drop=FALSE]
    ##res = res[grepl("GY:ID_",row.names(res)), , drop=FALSE]
    return(res)
}

getVariance = function(asr, comp){
	var = summary(asr)$varcomp
	idx = which(rownames(var)==comp)
	v = var$component[idx]
	print(paste("variance component", v))
	return(v)
}

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

# load data from centralized stations
load("data/diversity.panel.data.gp.rda")

###################
## get SNP data, convert it and impute it for use in genomic selection with rrBLUP

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

save(imputedcl, file="output/genoytpic.data.rrBLUP.rda")

#####################
## work on metric data

# check factor levels in metric data
head(met)
sapply(met, class)

f = c('LOCATION','PLOT','YEAR','REP','BLK','ROW','COL','ID')
head(met[f])

met[f] = lapply(met[f], as.factor)
sapply(met, class)

met[c('DB','DF','DM')] = lapply(met[c('DB','DF','DM')], as.numeric)
sapply(met, class)

head(met)

# produce BLUPs on metric data
trait = 'DB'
vioplot(met$DB ~ met$LOCATION + met$YEAR, col='lightgrey', ylab=trait)
colnames(met)
traits = colnames(met)[-c(1:8)] ## get list of traits

nLoc = 2
nRep = 2
nYear = 2

# set up dataframe for holding blups and h2
blups.met = data.frame(ID = paste("ID_",unique(met$ID), sep="") , row.names = unique(met$ID))
head(blups.met)
h2.met = data.frame()

# loop through traits with model for combined year, within year, within location, within location & year
for(t in traits){  ## t = traits[1]

  ## model for across year analysis  
  # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(LOCATION) + id(YEAR) + id(ID):id(LOCATION) + id(LOCATION):id(YEAR) + id(REP):id(LOCATION):id(YEAR), data=met, maxit=60, start.values=TRUE)
  # iv = asr$vparameters.table
  # iv$Constraint[]='U'
  # 
  asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
             random= ~ id(ID) 
             + id(LOCATION) 
             + id(YEAR) 
             + id(ID):id(LOCATION) 
             + id(ID):id(YEAR) 
             + id(LOCATION):id(YEAR) 
             + id(ID):id(LOCATION):id(YEAR) 
             + id(REP):id(LOCATION):id(YEAR)
             , data=met, maxit=100) ##, G.param = iv)
  
  print(summary(asr)$varcomp)
  
  
  b = getBLUPs(asr)
  b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", t)
  head(b.tmp)
      
  blups.met = merge(blups.met, b.tmp, by="ID")
  head(blups.met)

  
  
  Vg = getVariance(asr, 'ID')
  Vgl = getVariance(asr, 'ID:LOCATION')
  Vgy = getVariance(asr, 'ID:YEAR')
  Vgyl = getVariance(asr, 'ID:LOCATION:YEAR')
  Ve = getVariance(asr, 'units!R')
  
  h2 = Vg/(Vg + Vgl/nLoc + Vgy/nYear + Vgyl/(nLoc*nYear) + Ve/(nRep*nLoc*nYear)) ## h2 calculation for multiple locations in multiple years
  print(paste("H2 for", t, ":", h2))
  
  h2.met = rbind(h2.met, data.frame(trait=t, year="ALL", location="ALL", h2=h2))
  
  asr=NA
  
  ## models within year, across locations ##
  for(y in levels(met$YEAR)){  ## y=2012
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', y))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(ID):id(LOCATION) + id(LOCATION) + id(REP):id(LOCATION), data=met, subset= YEAR==y, maxit=20, start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'

    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(LOCATION) 
               + id(ID):id(LOCATION) 
               + id(REP):id(LOCATION), 
               data=met, subset= YEAR==y, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    name = paste(t, y, sep=".")
    b = getBLUPs(asr)
    b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
    head(b.tmp)
        
    blups.met = merge(blups.met, b.tmp, by="ID")
    head(blups.met)

    Vg = getVariance(asr, 'ID')
    Vgl = getVariance(asr, 'ID:LOCATION')
    Ve = getVariance(asr, 'units!R')

    h2 = Vg/(Vg + Vgl/nLoc + Ve/(nRep*nLoc))  ## calculation of H2, with one year at two or more locations
    print(paste("H2 for", t, "in year", y, ":", h2))
    
    h2.met = rbind(h2.met, data.frame(trait=t, year=y, location="ALL", h2=h2))

    asr=NA
  }   
  
  ## models for one location across years ##
    for(l in levels(met$LOCATION)){  ## l='geregera'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', l))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(ID):id(LOCATION) + id(LOCATION) + id(REP):id(LOCATION), data=met, subset= YEAR==y, maxit=20, start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'

    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(YEAR) 
               + id(ID):id(YEAR) 
               + id(REP):id(YEAR), 
               data=met, subset= LOCATION==l, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    name = paste(t, l, sep=".")
    b = getBLUPs(asr)
    b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
    head(b.tmp)
        
    blups.met = merge(blups.met, b.tmp, by="ID")
    head(blups.met)

    Vg = getVariance(asr, 'ID')
    Vgy = getVariance(asr, 'ID:YEAR')
    Ve = getVariance(asr, 'units!R')

    h2 = Vg/(Vg + Vgy/nYear + Ve/nRep) ## H2 calculation for one location in two or more years 
    print(paste("H2 for", t, "in", l, ":", h2))
    
    h2.met = rbind(h2.met, data.frame(trait=t, year="ALL", location=l, h2=h2))
    
    asr=NA
  }   
  
  ## models for single location-year ##
  for(y in levels(met$YEAR)){  ## t='DB'; y=2012; l='geregera'
    for(l in levels(met$LOCATION)){

      writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', y, 'in', l))
      
      asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
                 random= ~ id(ID) 
                 + id(REP), 
                 data=met, subset=c(YEAR==y & LOCATION==l), maxit=20)
     
       print(summary(asr)$varcomp)

      name = paste(t, y, l, sep=".")
      b = getBLUPs(asr)
      b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
      head(b.tmp)
          
      blups.met = merge(blups.met, b.tmp, by="ID")
      head(blups.met)
        
      
      Vg = getVariance(asr, 'ID')
      Ve = getVariance(asr, 'units!R')
      
      h2 = Vg/(Vg + Ve/(nRep))  ## H2 calculation for one location in one year
      print(paste("H2 for", t, "in", l, "in year", y, ":", h2))
      
      h2.met = rbind(h2.met, data.frame(trait=t, year=y, location=l, h2=h2))
      
      asr=NA
      
    }
  } 
}

round(h2.met$h2,2)
corrplot(cor(blups.met[,-1], use='complete.obs'))
corrplot(cor(blups.met[traits], use='complete.obs'))

# ## for standardizing to range of 0-10
#     d=''; if(mean(met[,t], na.rm=TRUE)>10){d='/10'}
#     
#     asr=asreml(fixed = as.formula(paste(t,d,'~ 1')), random= ~ id(ID) + id(ID):id(LOCATION) + id(LOCATION) + id(REP):id(LOCATION), data=met, subset= YEAR==y, maxit=20, start.values=TRUE)
#     

#####################
## work on farmer data

## collapse farmer data to key-value pairs
head(farm)
colnames(farm)[1]="LOCATION"  ## match to met data
f.names = colnames(farm)[-c(1:8)] ## get names of farmer traits

## transform data to format for mixed model ##
farm2 = gather(farm, key="farmer", value = 'SCORE', f.names)
head(farm2)

farm2 = separate(farm2, col='farmer', into=c('GROUP', 'TRAIT', 'FARMER'), sep="_")
head(farm2)

g = data.frame(GENDER = substring(farm2$GROUP,1,1)) ## add variable for gender
farm2 = cbind(farm2, g)

farm2 = spread(farm2, TRAIT, SCORE)
head(farm2)  

sapply(farm2, class) ## check if factors

n = colnames(farm2)[1:11]
farm2[n] = lapply(farm2[n], as.factor) ## change to factors
sapply(farm2, class)

traits = colnames(farm2)[-(1:11)]
farm2[traits] = lapply(farm2[traits], as.numeric) ## change phenotypes to numeric
sapply(farm2, class)

head(farm2)

## test of full model ## 

t = 'OVERALL'
asr=asreml(fixed = as.formula(paste(t,'~ 1')),
           random= ~ id(ID) 
           + id(LOCATION) 
           + id(ID):id(LOCATION) 
           + id(REP):id(LOCATION) 
           + id(GENDER) 
           + id(GENDER):id(LOCATION) 
           + id(GROUP):id(GENDER):id(LOCATION) 
           + id(FARMER):id(GROUP):id(GENDER):id(LOCATION) 
           + id(ID):id(GENDER) 
           + id(GENDER):id(LOCATION) 
           + id(ID):id(GENDER):id(LOCATION) 
           + id(ROW):id(LOCATION) 
           + id(COL):id(LOCATION) 
           + id(ID):id(ROW):id(COL):id(LOCATION),
           data=farm2, maxit=60)
  
print(summary(asr)$varcomp)

getBLUPs(asr)


traits = colnames(farm2)[-c(1:11)]

## set up dataframe for holding blups and h2
blups.farm = data.frame(ID = paste("ID_",unique(met$ID), sep="") , row.names = unique(met$ID))
head(blups.farm)

h2.farm = data.frame()

## loop through all traits

for(t in traits){  ## t = traits[1]

  ## model for across year analysis  
  # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID)  + id(LOCATION) + id(ID):id(LOCATION) + id(REP):id(LOCATION) + id(GENDER) + id(GENDER):id(LOCATION) + id(GROUP):id(GENDER):id(LOCATION) + id(FARMER):id(GROUP):id(GENDER):id(LOCATION) + id(ID):id(GENDER) + id(GENDER):id(LOCATION) + id(ID):id(GENDER):id(LOCATION) + id(ROW):id(LOCATION) + id(COL):id(LOCATION) + id(ID):id(ROW):id(COL):id(LOCATION), data=farm2, maxit=60, start.values=TRUE)
  # iv = asr$vparameters.table
  # iv$Constraint[]='U'
  # 
  asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
             random= ~ id(ID) 
             + id(LOCATION) 
             + id(ID):id(LOCATION) 
             + id(REP):id(LOCATION) 
             + id(GENDER) 
             + id(GENDER):id(LOCATION) 
             + id(GROUP):id(GENDER):id(LOCATION) 
             + id(FARMER):id(GROUP):id(GENDER):id(LOCATION) 
             + id(ID):id(GENDER) 
             + id(GENDER):id(LOCATION) 
             + id(ID):id(GENDER):id(LOCATION) 
             + id(ID):id(GROUP):id(GENDER):id(LOCATION) 
             + id(ROW):id(REP):id(LOCATION) 
             + id(COL):id(REP):id(LOCATION) 
             + id(ID):id(ROW):id(COL):id(LOCATION), 
             data=farm2, maxit=60) ##G.param=iv)
  
  print(summary(asr)$varcomp)
  
  
  Vg = getVariance(asr, 'ID')
  Vgl = getVariance(asr, 'ID:LOCATION')
  Vgm = getVariance(asr, 'ID:GENDER')
  Vgml = getVariance(asr, 'ID:GENDER:LOCATION')
  Vw = getVariance(asr, 'ID:ROW:COL:LOCATION')
  Ve = getVariance(asr, 'units!R')
  
  ## DON'T KNOW WHAT DO DO WITH FARMER AND GROUP VARIANCE ! .... NEED TO CHECK THIS CALCULATION
  h2 = Vg/(Vg + Vgl/nLoc + Vgm/2 + Vgml/(nLoc*2) + Vw/(nRep*nLoc) + Ve/(nRep*nLoc*30))  
  print(paste("H2 for", t, ":", h2))
  
  pred = predict(asr, 'ID')
  a = pred$avsed
  
  h2.farm = rbind(h2.farm, data.frame(trait=t, gender='both', year="2012", location="ALL", h2=h2, avsed=a))
  
  
  # ## simple model ##
  # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID)  + id(LOCATION) + id(ID):id(LOCATION) + id(REP):id(LOCATION), data=farm2, maxit=60) ##, G.param = iv)
  # print(summary(asr)$varcomp)
  # 
  # Vg = getVariance(asr, 'ID')
  # Vgl = getVariance(asr, 'ID:LOCATION')
  # Ve = getVariance(asr, 'units!R')
  # 
  # h2 = Vg/(Vg + Vgl/nLoc + Ve/(nRep*nLoc))
  # print(paste("H2 for", t, ":", h2))
  # 
  # pred = predict(asr, 'ID')
  # pred$avsed
  # 
  # 
  
  ## treating gender as fixed effect ##
  
  #   asr=asreml(fixed = as.formula(paste(t,'~ 1 + GENDER')), 
  #            random= ~ id(ID) 
  #            + id(LOCATION) 
  #            + id(ID):id(LOCATION) 
  #            + id(REP):id(LOCATION) 
  #            + id(GENDER):id(LOCATION) 
  #            + id(GROUP):id(GENDER):id(LOCATION) 
  #            + id(FARMER):id(GROUP):id(GENDER):id(LOCATION) 
  #            + id(ID):id(GENDER) 
  #            + id(GENDER):id(LOCATION) 
  #            + id(ID):id(GENDER):id(LOCATION) 
  #            + id(ID):id(GROUP):id(GENDER):id(LOCATION) 
  #            + id(ROW):id(REP):id(LOCATION) 
  #            + id(COL):id(REP):id(LOCATION) 
  #            + id(ID):id(ROW):id(COL):id(LOCATION), 
  #            data=farm2, maxit=60) ##G.param=iv)
  # 
  # print(summary(asr)$varcomp)
  # 
  
  
  
  
  b = getBLUPs(asr)
  b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", t)
  head(b.tmp)
      
  blups.farm = merge(blups.farm, b.tmp, by="ID")
  head(blups.farm)
  
  asr=NA


  ## BLUPs and H2 calculation for each location ##
  for(l in levels(farm2$LOCATION)){  ## l='hagreselam'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', l))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'

    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(REP) 
               + id(GENDER) 
               + id(GROUP):id(GENDER) 
               + id(ID):id(GENDER)
               + id(ID):id(GROUP):id(GENDER)
               + id(FARMER):id(GROUP):id(GENDER) 
               + id(ROW):id(REP) 
               + id(COL):id(REP)
               + id(ID):id(ROW):id(COL),
               data=farm2, subset= LOCATION==l, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    

    Vg = getVariance(asr, 'ID')
    Vgm = getVariance(asr, 'ID:GENDER')
    Vw = getVariance(asr, 'ID:ROW:COL')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgm/2 + Vw/(nRep) + Ve/(nRep*30)) ## NOT SURE, NEED TO CHECK THIS CALCULATION !!! ## 2 genders and 30 farmers
    print(paste("H2 for", t, 'in', l, ":", h2))
    
    pred = predict(asr, 'ID')
    a = pred$avsed
    
    h2.farm = rbind(h2.farm, data.frame(trait=t, gender='both', year='2012', location=l, h2=h2, avsed=a))
    
    
    
    # ## simple model ##
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
    # print(summary(asr)$varcomp)
    # Vg = getVariance(asr, 'ID')
    # Ve = getVariance(asr, 'units!R')
    # 
    # h2 = Vg/(Vg + Ve/nRep)  
    # print(paste("H2 for", t, "in year", y, ":", h2))
    # 
    # pred = predict(asr, 'ID')
    # pred$avsed
    
    ## blups 
    name = paste(t, l, sep=".")
    b = getBLUPs(asr)
    b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
    head(b.tmp)
        
    blups.farm = merge(blups.farm, b.tmp, by="ID")
    head(blups.farm)
    
    asr=NA

  }   
  

  ## BLUPs and H2 calculation for each gender ##
  for(g in levels(farm2$GENDER)){  ## l='hagreselam'; g='M'
    
    writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g))
    
    # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
    # iv = asr$vparameters.table
    # iv$Constraint[]='U'

    asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
               random= ~ id(ID) 
               + id(LOCATION) 
               + id(ID):id(LOCATION) 
               + id(REP):id(LOCATION) 
               + id(GROUP):id(LOCATION) 
               + id(FARMER):id(GROUP):id(LOCATION) 
               + id(ID):id(GROUP):id(LOCATION)
               + id(ROW):id(REP):id(LOCATION) 
               + id(COL):id(REP):id(LOCATION) 
               + id(ID):id(ROW):id(COL):id(LOCATION),
               data=farm2, subset=GENDER==g, maxit=60) ##, G.param = iv)
    
    print(summary(asr)$varcomp)
    
    Vg = getVariance(asr, 'ID')
    Vgl = getVariance(asr, 'ID:LOCATION')
    Vw = getVariance(asr, 'ID:ROW:COL:LOCATION')
    Ve = getVariance(asr, 'units!R')
    
    h2 = Vg/(Vg + Vgl/nLoc + Vw/(nRep*nLoc) + Ve/(nLoc*nRep*15))  ## 15 farmers for each gender
    print(paste("H2 for", t, 'for', g, ":", h2))
    
    pred = predict(asr, 'ID')
    a = pred$avsed

    
    h2.farm = rbind(h2.farm, data.frame(trait=t, gender=g, year='2012', location='ALL', h2=h2, avsed=a))
    

    ## blups
    name = paste(t, g, sep=".")
    b = getBLUPs(asr)
    b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
    head(b.tmp)
        
    blups.farm = merge(blups.farm, b.tmp, by="ID")
    head(blups.farm)
    
    asr=NA
    

  }     

  
  ## BLUPs and H2 calculation for each gender at each location ##
    for(l in levels(farm2$LOCATION)){  ## l='hagreselam', g='M'
      for(g in levels(farm2$GENDER)){
        
        writeLines("\n\n"); print(paste("RUNNING TRAIT:", t, 'for', g, 'in', l))
        
        # asr=asreml(fixed = as.formula(paste(t,'~ 1')), random= ~ id(ID) + id(REP) , data=farm2, subset= LOCATION==l, maxit=60) ##  start.values=TRUE)
        # iv = asr$vparameters.table
        # iv$Constraint[]='U'
        
        asr=asreml(fixed = as.formula(paste(t,'~ 1')), 
                   random= ~ id(ID) 
                   + id(REP) 
                   + id(GROUP) 
                   + id(FARMER):id(GROUP) 
                   + id(ID):id(GROUP)
                   + id(ROW):id(REP)
                   + id(COL):id(REP) 
                   + id(ID):id(ROW):id(COL),
                   data=farm2, subset= GENDER==g & LOCATION==l, maxit=60) ##, G.param = iv)
        
        print(summary(asr)$varcomp)
        
        
        Vg = getVariance(asr, 'ID')
        Vw = getVariance(asr, 'ID:ROW:COL')
        Ve = getVariance(asr, 'units!R')
        
        h2 = Vg/(Vg + Vw/(nRep) + Ve/(nRep*15))  ## 15 farmers for each gender
        print(paste("H2 for", t, 'for', g, 'in', l,  ":", h2))

        pred = predict(asr, 'ID')
        a = pred$avsed
        
        h2.farm = rbind(h2.farm, data.frame(trait=t, gender=g, year='2012', location=l, h2=h2, avsed=a))
        
        
        ## blups
        name = paste(t, l, g, sep=".")
        b = getBLUPs(asr)
        b.tmp = data.frame(row.names(b), b$effect); names(b.tmp) = c("ID", name)
        head(b.tmp)
        
        blups.farm = merge(blups.farm, b.tmp, by="ID")
        head(blups.farm)
        
        asr=NA
      }
  }     
}  

head(blups.farm)
h2.farm

#################
## save datasets

# Save the derived objects
save(blups.met, h2.met, blups.farm, h2.farm, file = "output/diversity.panel.BLUPs.rda")
