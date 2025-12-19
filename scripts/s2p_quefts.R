#' The QUEFTS model
#'
#' This function contains the QUEFTS model for NPK
#'
#' @param nin (numeric) the effective nitrogen input (kg N/ha) from soil and fertilizer
#' @param pin (numeric) the effective phosphorus input (kg P/ha) from soil and fertilizer
#' @param kin (numeric) the effective potassium input (kg K/ha) from soil and fertilizer
#' @param soc (numeric) the soil organic carbon content (g / kg)
#' @param pol (numeric) the Olsen-P content (mg P/kg)
#' @param phw (numeric) the soil pH measured in water
#' @param kexch (numeric) the exchangeable potassium content (mmol / kg)
#' @param country (character) The country for which to simulate, options: Kenya, Madagascar, Malawi, Mozambique, Tanzania, Uganda, Zambia
#' @param crop (character) the crop for which to simulate. Options: Maize, Potato
#' @param management_level (character) The level of intensity at which the field is managed, options: High, Medium, Low
#' @param max_yield (integer) the maximal attainable yield (kg / ha)
#' @param qwf (numeric) fraction affecting dilution vs accumulation
#'
#'
#'@export
s2p_calc_quefts <- function(nin,pin,kin, country, crop, management_level, max_yield, qwf = 2){
  
  # Visual bindings
  Country = Crop = Management = a1 = n1 = d1 = r1 = s1 = a2 = n2 = d2 = r2 = s2 = id = uptake = upt1 = upt2 = i.upt1 = YL = YD = YA = YM = yield = . = NULL
  
  # check input lenght
  arg.length <- max(length(nin),length(pin),length(kin),length(soc),length(pol),length(kexch),length(phw))
  
  # check input
  checkmate::check_numeric(nin,lower=0,upper = 500,len = arg.length)
  checkmate::check_numeric(pin,lower=0,upper = 500,len = arg.length)
  checkmate::check_numeric(kin,lower=0,upper = 800,len = arg.length)
  checkmate::check_numeric(pol,lower=0,upper = 100,len = arg.length)
  checkmate::check_numeric(soc,lower=0,upper = 200,len = arg.length)
  checkmate::check_numeric(phw,lower=3,upper = 8,len = arg.length)
  checkmate::check_numeric(kexch,lower=0,upper = 150,len = arg.length)
  
  # which QUEFTS are used as input
  nuts=c('N','P','K')
  
  # all combinations
  o1 <- utils::combn(nuts, 2)
  o1 <- as.data.table(t(cbind(o1, o1[c(2,1),])), stringsAsFactors = FALSE)
  colnames(o1) <- c('n1', 'n2')
  
  # Select quefts parameters for crop and management level
  parms_quefts <- S2P::parms_quefts[Country == country & Crop == crop & Management == management_level]
  
  # add QUEFTS properties to data.table
  o1[, a1 := as.numeric(parms_quefts[,mget(paste0('a',n1))])]
  o1[, d1 := as.numeric(parms_quefts[,mget(paste0('d',n1))])]
  o1[, r1 := as.numeric(parms_quefts[,mget(paste0('r',n1))])]
  o1[, a2 := as.numeric(parms_quefts[,mget(paste0('a',n2))])]
  o1[, d2 := as.numeric(parms_quefts[,mget(paste0('d',n2))])]
  o1[, r2 := as.numeric(parms_quefts[,mget(paste0('r',n2))])]
  
  # estimate soil supply of nutrients
  dsoil <- data.table(corg = corg,pols = pols,kexch = kexch, phw = phw)
  
    # estimate pH correction factors, so that within 0-1 range (Janssen, 1998)
    dsoil[,fphn := pmin(1,pmax(0,0.25 * (phw - 3)))]
    dsoil[,fphp := pmin(1,pmax(0,1 - 0.5 * (phw-6)^2))]
    dsoil[,fphk := pmin(1,pmax(0,0.625 * (3.4 - 0.4 * phw)))]
    
    # estimate the soil nutrient supply
    dsoil[]
    
  # set internal data.table with the NPK input (kg / ha)
  dt0 <- data.table(N = nin,P = pin,K = kin)
  
  # set internal table with maximal yields per field
  dt.ym <- data.table(YM = max_yield)
  
  # set an unique ID per combination of inputs
  dt0[,id:=.I]
  dt.ym[,id := .I]
  
  # spread the data.table with input data for joining purposes later
  dt <- melt(dt0,
             id.vars='id',
             measure.vars=nuts,
             variable.name='n1',
             value.name='s1')
  
  # enumerate quefts parameter set with all 2-way interactions for each grid cel
  o2 = o1[rep(1:nrow(o1),nrow(dt)),]
  o2[,id := sort(rep(min(dt$id):max(dt$id),3 * nrow(o1)))]
  setkey(o2,id,n1,n2)
  setkey(dt,id,n1)
  
  # left_join with element n1 and n2
  o2[dt,on=.(id,n1),s2:=s1]
  o2=o2[dt,on=c(id='id',n2='n1'),nomatch=0L]
  setnames(o2,names(o2),c(names(o2)[1:9],'s1','s2'))
  
  # calculate for each soil the uptake for combinations of nutrients, and then the minimum of two-pair interactions
  o2[, uptake := qnup(a1,d1,r1,s1+0.001,a2,d2,r2,s2+0.001)][,upt1 := min(uptake),by=list(id,n1)]
  
  # add minimum yield corresponding to second nutrient (necessary for calculation YA and YD later)
  o2[o2[,c('id','n1','upt1')],on=c(id='id',n2='n1'),upt2 := i.upt1]
  up = dcast(o2,id~n1,value.var='upt1',fun.aggregate=mean)
  
  # calculate yields at accumulation (YA.) and dilution (YD.) per nutrient (DELAYING)
  o2[,c('YA','YD','YL','a1','d1','r1'):= qYAYD(n1,n2,upt2,upt1,qp=as.data.frame(parms_quefts))]
  
  # merge table with maximal yields per field
  o2 <- merge(o2,dt.ym, by = 'id')
  
  # calculate maximum yield, and yield per pair of nutrients
  o2[,c('YL'):= list(round(min(YL,YM))),by=id][,c('yield'):=list(qYpNpair(upt1,YD,YA,a1,d1,r1,YL))]
  
  # y1 equals the average of all pair-estimates
  y1 <- o2[,list(yield=mean(yield,na.rm=TRUE)),by=id]
  
  # merge with the estimated nutrient uptake
  y1 <- y1[up,on='id']
  setnames(y1,names(y1),ifelse(names(y1) %in% nuts,paste0('u',names(y1)),names(y1)))
  
  # ratio yield vs uptake
  score1 = as.matrix(y1$yield / y1[,mget(paste0('u',nuts))],ncol=length(nuts))
  # subtract by average dilution and accumulation, (a+d)/2
  score2 = sweep(score1,2,as.matrix(parms_quefts[,mget(paste0('C',nuts))]/qwf))
  # make score 2 relative
  score3 = sweep(score2,2,FUN ='/',as.matrix(parms_quefts[,mget(paste0('C',nuts))]/qwf))
  # sum error for all nutrients, and weigh the errors
  score4 = score3^2 %*% t(parms_quefts[,mget(paste0('s',nuts))])
  # add score to the output
  y1$score = c(score4)
  
  # add the given input dose
  y1 = y1[dt0,on='id']
  
  return(y1)
}



#' Derive nutrient requirement to obtain target yield and optimum nutrient content in crop
#'
#' This function predicts the uptake and yield for a given simulation output object
#'
#' @param inp (data.table)
#' @param inp2 (integer)
#' @param nuts (integer)
#' @param qdb1 (integer)
#' @param qyield (integer)
#'
#' @details
#' predict the uptake and yield for a given simulation output object (unfertilized and fertilized)
#'
#' @import data.table
#'
#' @export
qFindNutReq <- function(inp,inp2,nuts,qdb1,qyield){

  # Visual bindings
  . = tyield_scaled = runid = scrore = field_id = gId = type = score = NULL

  # situations where crop yield in unfertilized situation is higher than required, then requirement = soil supply
  inp0s <- inp2[tyield_scaled<=qyield$ini_yield,]
  inp0s <- inp0s[qyield,on=c('field_id'='field_id')][!is.na(tyield_scaled)]
  inp0s <- inp0s[,mget(c('field_id','tyield_scaled',paste0('av',rep(nuts,2)),paste0('ini_',nuts),'ini_score'))]
  setnames(inp0s,c('field_id','tyield_scaled',paste0('av',nuts),nuts,paste0('u',nuts),'score'))

  # local copy of inp2 en qdb2
  inp2s <- inp2[tyield_scaled>qyield$ini_yield,]
  qdb2  <- qdb1[tyield_scaled %in% inp2$tyield_scaled, ]#[,runid := NULL]

  # setkey sorts the table in ascending order
  setkey(inp2s,tyield_scaled)
  setkeyv(qdb2,c('tyield_scaled','score',nuts))

  # conditional join minimum nutrient supply needed for target yield where dose is higher than soil supply
  # inp2s[, c(nuts, 'score', paste0('u', nuts)) :=
  #         # perform the join, and retrieve required fertilizer dose (and score)
  #         merge(inp2s, qdb2[, c(paste0('f', nuts), 'score', paste0('u', nuts)),
  #                           on = c('tyield_scaled', paste0(nuts, ">=av", nuts)),
  #                           mult = 'first'],
  #               by = c(nuts),
  #               suffixes = c('', '.qdb'))[, ..c(nuts, 'score', paste0('u', nuts))]
  # ]


  inp2s[, c(nuts,'score',paste0('u',nuts)):=
          # perform the join, and retrieve required fertilizer dose (and score)
          qdb2[inp2s, mget(c(paste0('f',nuts),'score',paste0('u',nuts))),
               # join conditions, and get the first row meeting the conditions (and given setkey order)
               on=c('tyield_scaled',paste0(nuts,">=av",nuts)),mult='first']
  ]

  # remove rows that have no score or id
  inp2s <- inp2s[!is.na(score),][!is.na(field_id),]

  # retreive minimum nutrient supply for all nutrients (kg/ ha) to obtain target yields (subset of qdb1)
  qdb3 <- qdb2[qdb2[,.I[which.min(score)],by=tyield_scaled]$V1]

  # select samples where nutrient supply is higher than required, and add minimum dose required to obtain yield with optimum nutrient content
  inp3s <- inp2[!field_id %in% c(inp0s$field_id,inp2s$field_id),mget(c('field_id','tyield_scaled',paste0('av',nuts)))][qdb3[,mget(c(nuts,'score','tyield_scaled',paste0('u',nuts)))],on=.(tyield_scaled)]

  # add default filter for all elements (nutrient supply should be higher than soil supply)
  inp3s[,c(paste0('sf',nuts)):=lapply(nuts,function(i) paste0(i,'>=av',i))]

  # change the default filter for relevant cases
  for(i in nuts){inp3s[eval(parse(text=paste0(i,'<','av',i))),paste0('sf',i):=paste0(i,'<av',i)]}

  # add for each combined filter group (for all nutrients) an unique grouping number
  inp3s[,gId := .GRP,by=c(paste0('sf',nuts))]

  # which are the unique filters
  filt_unique <- unique(inp3s[,mget(paste0('sf',nuts))])

  # select relevant columns
  inp3s <- inp3s[,mget(c('field_id','tyield_scaled',paste0('av',nuts),'gId'))]

  # add nutrient supply for those grid cells where one of the soils' supply is higher than required
  inp3s[,c(nuts,'score',paste0('u',nuts)):= qYsel(sgId=gId,sdb=inp3s,snuts=nuts,sqdb=qdb2,sfilt=filt_unique),by=gId]
  inp3s <- inp3s[!is.na(field_id),]

  # replace nutrient requirement with the maximum of soil supply and nutrient needed
  for(i in nuts){inp3s[,paste0(i):= do.call(pmax,.SD),.SDcols=c(i,paste0('av',i))]}

  # check for id's that have no score and add minimum
  inp4s <- inp2[!field_id %in% c(inp0s$field_id,inp2s$field_id,inp3s$field_id),mget(c('field_id','tyield_scaled',paste0('av',nuts)))][qdb3[,mget(c(nuts,'score','tyield_scaled',paste0('u',nuts)))],on=.(tyield_scaled)]
  inp4s <- inp4s[!is.na(field_id)]

  # in the case that the expected yield is lower than minimum of quefts, replace with minumum of quefts
  inp5s <- inp2[!field_id %in% c(inp0s$field_id,inp2s$field_id,inp3s$field_id,inp4s$field_id),mget(c('field_id','tyield_scaled'))]
  inp5s[,tyield_scaled := pmax(tyield_scaled,min(qdb1$tyield_scaled))]
  inp5s <- inp5s[qdb3[,mget(c(nuts,'score','tyield_scaled',paste0('u',nuts)))],on=.(tyield_scaled)]
  inp5s <- inp5s[!is.na(field_id)]

  # add type recommendation (nutrient requirement is always higher than soil supply (1) otherwise (2))
  inp0s[,type:=0]
  inp2s[,type:=1]
  inp3s[,type:=2]
  inp4s[,type:=3]
  inp5s[,type:=4]

  # combine the two databases and combine with original soil properties
  svars = c('field_id','tyield_scaled',nuts,'score',paste0('u',nuts),'type')
  #qm <- if(nrow(inp3s)>0){rbind(inp2s[,mget(svars)],inp3[,mget(svars)])} else {inp2s[,mget(svars)]}
  qm <- rbind(inp0s[,mget(svars)],inp2s[,mget(svars)],inp3s[,mget(svars)],inp4s[,mget(svars)],inp5s[,mget(svars)])
  setnames(qm,nuts,paste0('req',nuts))

  # add select keys to order both data.tables,
  setkey(qm,field_id,tyield_scaled)
  setkey(inp,field_id,tyield_scaled)
  qm <- inp[qm,on=.(field_id,tyield_scaled)]

  return(qm)
  }
