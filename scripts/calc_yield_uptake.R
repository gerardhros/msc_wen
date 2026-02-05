# calculate yield from Nutrient uptake

# require package
require(data.table)

# which QUEFTS are used as input
nuts=c('N','P','K')

# make a table with all combinations of selected nutrients
o1 <- utils::combn(nuts, 2)
o1 <- as.data.table(t(cbind(o1, o1[c(2,1),])), stringsAsFactors = FALSE)
colnames(o1) <- c('n1', 'n2')

# Select quefts parameters for crop and management level
parms_quefts <- fread('C:/ESA/02 msc projects/rui wen/msc_wen/data/parms_quefts.csv')
parms_quefts <- parms_quefts[Country == 'Uganda' & Crop == 'Maize' & Management == 'High']

# add QUEFTS properties to data.table with all nutrient combinations
o1[, a1 := as.numeric(parms_quefts[,mget(paste0('a',n1))])]
o1[, d1 := as.numeric(parms_quefts[,mget(paste0('d',n1))])]
o1[, r1 := as.numeric(parms_quefts[,mget(paste0('r',n1))])]
o1[, a2 := as.numeric(parms_quefts[,mget(paste0('a',n2))])]
o1[, d2 := as.numeric(parms_quefts[,mget(paste0('d',n2))])]
o1[, r2 := as.numeric(parms_quefts[,mget(paste0('r',n2))])]

# add the supply of soil plus nutrients (kg K/ha) for each combination of nutrients
o1[n1=='K',s1 := 175]
o1[n2=='K',s2 := 175]
o1[n1=='N',s1 := 220]
o1[n2=='N',s2 := 220]
o1[n1=='P',s1 := 54]
o1[n2=='P',s2 := 54]

# add id for the field (only single field here)
o1[,id:=1]

# calculate for each soil the uptake for combinations of nutrients, and then the minimum of two-pair interactions
uptake <- o1[, uptake := qnup(a1,d1,r1,s1+0.001,a2,d2,r2,s2+0.001)][,min(uptake),by=list(n1)]

  # save the uptake for later use
  up <- dcast(uptake,.~n1,value.var='V1')[,.:=NULL][,id:=1]
  
# add this for both first and second nutrient in the combination
o1 <- merge(o1,uptake[,.(n1,upt1 = V1)],by='n1')
o1 <- merge(o1,uptake[,.(n2 = n1,upt2 = V1)],by='n2')

# calculate yields at accumulation (YA.) and dilution (YD.) per nutrient (DELAYING)
o1[,c('YA','YD','YL','a1','d1','r1'):= qYAYD(n1,n2,upt2,upt1,qp=as.data.frame(parms_quefts))]

# add a maximum yield
o1[, YM := 30000]

# calculate maximum yield, and yield per pair of nutrients
o1[,'YL':= list(round(min(YL,YM))),by=id]
o1[,'yield':=list(qYpNpair(upt1,YD,YA,a1,d1,r1,YL))]

# calculate yield equalling the average of all pair-estimates
y1 <- o1[,list(yield=mean(yield,na.rm=TRUE)),by=id]

# merge with the estimated nutrient uptake
y1 <- merge(y1,up,by='id')
setnames(y1,names(y1),ifelse(names(y1) %in% nuts,paste0('u',names(y1)),names(y1)))

# function to calculate the yield as funciton of uptake
calc_yield_from_uptake <- function(nup,kup,pup){
  
  # make internal table for a single site (id = 1)
  d1 <- data.table(N = nup, P = pup, K = kup, id = 1)
  
  # melt the datatable
  uptake <- melt(d1,id.vars='id',variable.name = 'n1',value.name = 'upt')
  uptake[,n1 := as.character(n1)]
  
  # which QUEFTS are used as input
  nuts=c('N','P','K')
  
  # make a table with all combinations of selected nutrients
  o1 <- utils::combn(nuts, 2)
  o1 <- as.data.table(t(cbind(o1, o1[c(2,1),])), stringsAsFactors = FALSE)
  colnames(o1) <- c('n1', 'n2')
  
  # add this for both first and second nutrient in the combination
  o1 <- merge(o1,uptake[,.(id,n1,upt1 = upt)],by='n1')
  o1 <- merge(o1,uptake[,.(n2 = n1,upt2 = upt)],by='n2')
  
  # calculate yields at accumulation (YA.) and dilution (YD.) per nutrient (DELAYING)
  o1[,c('YA','YD','YL','a1','d1','r1'):= qYAYD(n1,n2,upt2,upt1,qp=as.data.frame(parms_quefts))]
  
  # add a maximum yield
  o1[, YM := 30000]
  
  # calculate maximum yield, and yield per pair of nutrients
  o1[,'YL':= list(round(min(YL,YM))),by=id]
  o1[,'yield':=list(qYpNpair(upt1,YD,YA,a1,d1,r1,YL))]
  
  # calculate yield equalling the average of all pair-estimates
  y1 <- o1[,list(yield=mean(yield,na.rm=TRUE)),by=id]
  
  # return yield as number
  y1 <- y1[,yield]
  
  return(y1)
}
