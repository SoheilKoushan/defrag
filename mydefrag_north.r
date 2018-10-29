#
# Trim RDS files, add key columns, and defragment
#
#
# Plot segim for individual object
#
#####################
# library on PAWSEY #
#####################
.libPaths(c('/group/pawsey0160/skoushan/R_Library/',.libPaths()))
library(EBImage)
library(celestial)
library(devtools)
#
.libPaths(c('/group/pawsey0160/skoushan/Zeus/Libs/',.libPaths()))
library(Cairo)
library(ProFound)
library('magicaxis')
library('data.table')
library(dplyr)
require(foreign)
require(MASS)
##################
# Reference Data #
##################
allz=fread("/group/pawsey0160/waves/ref/AllGAMACounts/allzcounts.csv")
ebv=fread("/group/pawsey0160/waves/ref/ebvtrim.csv")
gaia = fread('/group/pawsey0160/waves/ref/gaiastarmaskwaves.csv')
gamaprior=fread('/group/pawsey0160/waves/ref/Specmag.csv')
gaiaraoff=0.0
gaiadecoff=0.0
############################### 
# Read in ProFound RDS output #
###############################
for (i in seq(167.5,169.5)){ # seq(128.5,237.5)     
  for (j in seq(-3.5,3.5)){     

		Dir = '/group/pawsey0160/waves/rds/'
		if (file.exists(paste(Dir,'RA',i,'DEC',j,'/','waves_pro_stack.rds',sep='')) == TRUE){
		everything=readRDS(paste(Dir,'RA',i,'DEC',j,'/','waves_pro_stack.rds',sep=''))
		datafile <- as.data.table(cbind(everything$pro_detect$segstats,everything$cat_col,everything$cat_tot))
        } 

# SOUTH
#for (i in seq(-31.5,53.5)){
 # for (j in seq(-35.5,-26.5)){

  #		ifelse(i<0,i%%360,i) #if (i < 0 ){i = i%%360}

	#	Dir = '/group/pawsey0160/waves/rds/'
	#	everything=readRDS(paste(Dir,'RA',i,'DEC',j,'/','waves_pro_stack.rds',sep=''))
	#	datafile <- as.data.table(cbind(everything$pro_detect$segstats,everything$cat_col,everything$cat_tot))

	 # }
	#}		        
##################################
# Make a list of ProFound output #
##################################
# determine whether a row exists within a data frame or not
if (nrow(merge(everything$cat_tot$mag_ut[1],everything$cat_tot)) != 0){
datafile0=datafile[,list(uniqueID,segID,xmax,ymax,xcen,ycen,RAcen,Deccen,RAmax,Decmax,flux_ut,flux_err_ut,flux_gt,flux_err_gt,flux_rt,flux_err_rt,flux_it,flux_err_it,flux_Zt,flux_err_Zt,flux_Yt,flux_err_Yt,flux_Jt,flux_err_Jt,flux_Ht,flux_err_Ht,flux_Kt,flux_err_Kt,flux_uc,flux_err_uc,flux_gc,flux_err_gc,flux_rc,flux_err_rc,flux_ic,flux_err_ic,flux_Zc,flux_err_Zc,flux_Yc,flux_err_Yc,flux_Jc,flux_err_Jc,flux_Hc,flux_err_Hc,flux_Kc,flux_err_Kc,sky_mean,skyRMS_mean,SB_N50_Zt,SB_N100_Zt,SB_N50_rt,SB_N100_rt,R50,R90,R100,N100,mag_uc,mag_ut,mag_gc,mag_gt,mag_rc,mag_rt,mag_ic,mag_it,mag_Zc,mag_Zt,mag_Yc,mag_Yt,mag_Jc,mag_Jt,mag_Hc,mag_Ht,mag_Kc,mag_Kt)]
} else {(datafile0=datafile[,list(uniqueID,segID,xmax,ymax,xcen,ycen,RAcen,Deccen,RAmax,Decmax,flux_Zt,flux_err_Zt,flux_Yt,flux_err_Yt,flux_Jt,flux_err_Jt,flux_Ht,flux_err_Ht,flux_Kt,flux_err_Kt,flux_Zc,flux_err_Zc,flux_Yc,flux_err_Yc,flux_Jc,flux_err_Jc,flux_Hc,flux_err_Hc,flux_Kc,flux_err_Kc,sky_mean,skyRMS_mean,SB_N50_Zt,SB_N100_Zt,R50,R90,R100,N100,mag_Zc,mag_Zt,mag_Yc,mag_Yt,mag_Jc,mag_Jt,mag_Hc,mag_Ht,mag_Kc,mag_Kt)])}
#####################################
# Add Name info and estimate seeing #
#####################################
 iauname=IAUID(datafile0$RAcen,datafile0$Deccen,name = 'WAVES',epoch = '_N_')
 frameid=ceiling(min(datafile0$RAcen))*100.0+ceiling(min(datafile0$Deccen))
 uberid=frameid*1e10+datafile0$uniqueID
 seeing=log10(median(datafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50]))
 datafile0=cbind(iauname,frameid,uberid,datafile0,seeing)
 names(datafile0)[1] <- "IAUID"
 names(datafile0)[2] <- "FrameID"
 names(datafile0)[3] <- "uberID"
 names(datafile0)[length(datafile0)] <- "seeing"
#################################################################
# Trim Dustmask and Starmask to appropriate region (saves time) #
#################################################################
wavesra=c(min(datafile0[,RAcen]),max(datafile0[,RAcen]))
wavesdec=c(min(datafile0[,Deccen]),max(datafile0[,Deccen]))
gaia0=gaia[ra > wavesra[1]-1.0 & ra < wavesra[2]+1.0 & dec > wavesdec[1]-1.0 & dec < wavesdec[2]+1.0,]
ebv0=ebv[RA_J2000 > wavesra[1]-1.0 & RA_J2000 < wavesra[2]+1.0 & DEC_J2000 > wavesdec[1]-1.0 & DEC_J2000 < wavesdec[2]+1.0,]
##################################################
# Determine extinction (E(B-V) using Planck map) #
##################################################
waves=coordmatch(coordref=datafile0[,list(RAmax,Decmax)],coordcompare=ebv0[,list(RA_J2000,DEC_J2000)],rad=150,inunitref="deg",inunitcompare="deg",radunit="asec")
datafile0=cbind(datafile0,ebv0[waves$ID[,1],"EBV"])
####################################################################################
# Correct all mags for attenuation using coefficients from Liske et al (2015)      #
#                                                                                  #
#                                                                                  #
# Correct all mags for attenuation using coefficients from Galactic Extinction DMU #
####################################################################################
optmult=10^{0.4*(8.9-0.0)}
nirmult=10^{0.4*(8.9-30.0)}
#  
Au=5.155
Ag=3.793
Ar=2.751
Ai=2.086
AZ=1.71
AY=1.31
AJ=0.928
AH=0.592
AK=0.388
#
datafile0$mag_ut=ifelse(datafile0$mag_ut==-999,-999,datafile0$mag_ut-Au*datafile0$EBV)
datafile0$mag_gt=ifelse(datafile0$mag_gt==-999,-999,datafile0$mag_gt-Ag*datafile0$EBV)
datafile0$mag_rt=ifelse(datafile0$mag_rt==-999,-999,datafile0$mag_rt-Ar*datafile0$EBV)
datafile0$mag_it=ifelse(datafile0$mag_it==-999,-999,datafile0$mag_it-Ai*datafile0$EBV)
datafile0$mag_Zt=ifelse(datafile0$mag_Zt==-999,-999,datafile0$mag_Zt-AZ*datafile0$EBV)
datafile0$mag_Yt=ifelse(datafile0$mag_Yt==-999,-999,datafile0$mag_Yt-AY*datafile0$EBV)
datafile0$mag_Jt=ifelse(datafile0$mag_Jt==-999,-999,datafile0$mag_Jt-AJ*datafile0$EBV)
datafile0$mag_Ht=ifelse(datafile0$mag_Ht==-999,-999,datafile0$mag_Ht-AH*datafile0$EBV)
datafile0$mag_Kt=ifelse(datafile0$mag_Kt==-999,-999,datafile0$mag_Kt-AK*datafile0$EBV)
#
datafile0$mag_uc=ifelse(datafile0$mag_uc==-999,-999,datafile0$mag_uc-Au*datafile0$EBV)
datafile0$mag_gc=ifelse(datafile0$mag_gc==-999,-999,datafile0$mag_gc-Ag*datafile0$EBV)
datafile0$mag_rc=ifelse(datafile0$mag_rc==-999,-999,datafile0$mag_rc-Ar*datafile0$EBV)
datafile0$mag_ic=ifelse(datafile0$mag_ic==-999,-999,datafile0$mag_ic-Ai*datafile0$EBV)
datafile0$mag_Zc=ifelse(datafile0$mag_Zc==-999,-999,datafile0$mag_Zc-AZ*datafile0$EBV)
datafile0$mag_Yc=ifelse(datafile0$mag_Yc==-999,-999,datafile0$mag_Yc-AY*datafile0$EBV)
datafile0$mag_Jc=ifelse(datafile0$mag_Jc==-999,-999,datafile0$mag_Jc-AJ*datafile0$EBV)
datafile0$mag_Hc=ifelse(datafile0$mag_Hc==-999,-999,datafile0$mag_Hc-AH*datafile0$EBV)
datafile0$mag_Kc=ifelse(datafile0$mag_Kc==-999,-999,datafile0$mag_Kc-AK*datafile0$EBV)
#
datafile0$flux_ut=ifelse(datafile0$flux_ut==0.0,-999,datafile0$flux_ut*10^(0.4*(Au*datafile0$EBV))*optmult)
datafile0$flux_gt=ifelse(datafile0$flux_gt==0.0,-999,datafile0$flux_gt*10^(0.4*(Ag*datafile0$EBV))*optmult)
datafile0$flux_rt=ifelse(datafile0$flux_rt==0.0,-999,datafile0$flux_rt*10^(0.4*(Ar*datafile0$EBV))*optmult)
datafile0$flux_it=ifelse(datafile0$flux_it==0.0,-999,datafile0$flux_it*10^(0.4*(Ai*datafile0$EBV))*optmult)
datafile0$flux_Zt=ifelse(datafile0$flux_Zt==0.0,-999,datafile0$flux_Zt*10^(0.4*(AZ*datafile0$EBV))*nirmult)
datafile0$flux_Yt=ifelse(datafile0$flux_Yt==0.0,-999,datafile0$flux_Yt*10^(0.4*(AY*datafile0$EBV))*nirmult)
datafile0$flux_Jt=ifelse(datafile0$flux_Jt==0.0,-999,datafile0$flux_Jt*10^(0.4*(AJ*datafile0$EBV))*nirmult)
datafile0$flux_Ht=ifelse(datafile0$flux_Ht==0.0,-999,datafile0$flux_Ht*10^(0.4*(AH*datafile0$EBV))*nirmult)
datafile0$flux_Kt=ifelse(datafile0$flux_Kt==0.0,-999,datafile0$flux_Kt*10^(0.4*(AK*datafile0$EBV))*nirmult)
#
datafile0$flux_err_ut=ifelse(datafile0$flux_err_ut==0.0,-999,datafile0$flux_err_ut*10^(0.4*(Au*datafile0$EBV))*optmult)
datafile0$flux_err_gt=ifelse(datafile0$flux_err_gt==0.0,-999,datafile0$flux_err_gt*10^(0.4*(Ag*datafile0$EBV))*optmult)
datafile0$flux_err_rt=ifelse(datafile0$flux_err_rt==0.0,-999,datafile0$flux_err_rt*10^(0.4*(Ar*datafile0$EBV))*optmult)
datafile0$flux_err_it=ifelse(datafile0$flux_err_it==0.0,-999,datafile0$flux_err_it*10^(0.4*(Ai*datafile0$EBV))*optmult)
datafile0$flux_err_Zt=ifelse(datafile0$flux_err_Zt==0.0,-999,datafile0$flux_err_Zt*10^(0.4*(AZ*datafile0$EBV))*nirmult)
datafile0$flux_err_Yt=ifelse(datafile0$flux_err_Yt==0.0,-999,datafile0$flux_err_Yt*10^(0.4*(AY*datafile0$EBV))*nirmult)
datafile0$flux_err_Jt=ifelse(datafile0$flux_err_Jt==0.0,-999,datafile0$flux_err_Jt*10^(0.4*(AJ*datafile0$EBV))*nirmult)
datafile0$flux_err_Ht=ifelse(datafile0$flux_err_Ht==0.0,-999,datafile0$flux_err_Ht*10^(0.4*(AH*datafile0$EBV))*nirmult)
datafile0$flux_err_Kt=ifelse(datafile0$flux_err_Kt==0.0,-999,datafile0$flux_err_Kt*10^(0.4*(AK*datafile0$EBV))*nirmult)
#
datafile0$flux_uc=ifelse(datafile0$flux_uc==0.0,-999,datafile0$flux_uc*10^(0.4*(Au*datafile0$EBV))*optmult)
datafile0$flux_gc=ifelse(datafile0$flux_gc==0.0,-999,datafile0$flux_gc*10^(0.4*(Ag*datafile0$EBV))*optmult)
datafile0$flux_rc=ifelse(datafile0$flux_rc==0.0,-999,datafile0$flux_rc*10^(0.4*(Ar*datafile0$EBV))*optmult)
datafile0$flux_ic=ifelse(datafile0$flux_ic==0.0,-999,datafile0$flux_ic*10^(0.4*(Ai*datafile0$EBV))*optmult)
datafile0$flux_Zc=ifelse(datafile0$flux_Zc==0.0,-999,datafile0$flux_Zc*10^(0.4*(AZ*datafile0$EBV))*nirmult)
datafile0$flux_Yc=ifelse(datafile0$flux_Yc==0.0,-999,datafile0$flux_Yc*10^(0.4*(AY*datafile0$EBV))*nirmult)
datafile0$flux_Jc=ifelse(datafile0$flux_Jc==0.0,-999,datafile0$flux_Jc*10^(0.4*(AJ*datafile0$EBV))*nirmult)
datafile0$flux_Hc=ifelse(datafile0$flux_Hc==0.0,-999,datafile0$flux_Hc*10^(0.4*(AH*datafile0$EBV))*nirmult)
datafile0$flux_Kc=ifelse(datafile0$flux_Kc==0.0,-999,datafile0$flux_Kc*10^(0.4*(AK*datafile0$EBV))*nirmult)
#
datafile0$flux_err_uc=ifelse(datafile0$flux_err_uc==0.0,-999,datafile0$flux_err_uc*10^(0.4*(Au*datafile0$EBV))*optmult)
datafile0$flux_err_gc=ifelse(datafile0$flux_err_gc==0.0,-999,datafile0$flux_err_gc*10^(0.4*(Ag*datafile0$EBV))*optmult)
datafile0$flux_err_rc=ifelse(datafile0$flux_err_rc==0.0,-999,datafile0$flux_err_rc*10^(0.4*(Ar*datafile0$EBV))*optmult)
datafile0$flux_err_ic=ifelse(datafile0$flux_err_ic==0.0,-999,datafile0$flux_err_ic*10^(0.4*(Ai*datafile0$EBV))*optmult)
datafile0$flux_err_Zc=ifelse(datafile0$flux_err_Zc==0.0,-999,datafile0$flux_err_Zc*10^(0.4*(AZ*datafile0$EBV))*nirmult)
datafile0$flux_err_Yc=ifelse(datafile0$flux_err_Yc==0.0,-999,datafile0$flux_err_Yc*10^(0.4*(AY*datafile0$EBV))*nirmult)
datafile0$flux_err_Jc=ifelse(datafile0$flux_err_Jc==0.0,-999,datafile0$flux_err_Jc*10^(0.4*(AJ*datafile0$EBV))*nirmult)
datafile0$flux_err_Hc=ifelse(datafile0$flux_err_Hc==0.0,-999,datafile0$flux_err_Hc*10^(0.4*(AH*datafile0$EBV))*nirmult)
datafile0$flux_err_Kc=ifelse(datafile0$flux_err_Kc==0.0,-999,datafile0$flux_err_Kc*10^(0.4*(AK*datafile0$EBV))*nirmult)
################################################
# Add blank classification columns to datafile #
################################################
wavesfile=data.frame()
datafile0 = as.data.table(cbind(datafile0,censep=as.numeric(sqrt((11000/2.0-datafile0[,"xmax"])^2+(11000/2.0-datafile0[,"ymax"])^2)),RAGAIA=datafile0$RAcen-gaiaraoff,DecGAIA=datafile0$Deccen-gaiadecoff,class="notclassified",duplicate=1,starscol=0,starssize=0,mask=0,noOPT=0,noIR=0,regroup=0,starmask=0,CATAID=0,Z=-999,NQ=0,PROB=0))
#
# Assign GAMA CATAIDS from input catalogues
#
#radius=1.5
#gamamatch=coordmatch(coordcompare=datafile0[,list(RAmax,Decmax)],coordref=gamainputcat[,list(RA,DEC)],rad=radius,inunitref="deg",inunitcompare="deg",radunit="asec")
#datafile0[gamamatch$bestmatch$compareID,"CATAID"]=gamainputcat[gamamatch$bestmatch$refID,"CATAID"]
##################################
# Assign regions with no IR data #
##################################
datafile0[flux_rt==-999,"noOPT"]=1
datafile0[flux_Zt==-999,"noIR"]=1
###################################
# Assign Stars by colour and size #
###################################
if (is.na(datafile0$mag_rt[1])==FALSE){
datafile0[,"starscol"]=0
datafile0[mag_rt > 19 & (mag_Jc-mag_Kc) < (0.05+0.05*(mag_rt-19)),"starscol"]=1
datafile0[mag_rt < 19 & (mag_Jc-mag_Kc) < 0.05 | mag_rt < 12,"starscol"]=3
datafile0[mag_rt > 19 & (mag_Jc-mag_Kc) < (0.05-0.1*(mag_rt-19)^2.0),"starscol"]=3
#
datafile0[,"starssize"]=0
datafile0[log10(R50) < seeing+0.1,"starssize"]=1
datafile0[log10(R50) < (seeing+0.1-0.15*(mag_rt-19.0)^1.0),"starssize"]=3
} else {
	datafile0[,"starscol"]=0
	datafile0[mag_Zt > 19 & (mag_Jc-mag_Kc) < (0.05+0.05*(mag_Zt-19)),"starscol"]=1
	datafile0[mag_Zt < 19 & (mag_Jc-mag_Kc) < 0.05 | mag_Zt < 12,"starscol"]=3
	datafile0[mag_Zt > 19 & (mag_Jc-mag_Kc) < (0.05-0.1*(mag_Zt-19)^2.0),"starscol"]=3
	#
	datafile0[,"starssize"]=0
	datafile0[log10(R50) < seeing+0.1,"starssize"]=1
	datafile0[log10(R50) < (seeing+0.1-0.15*(mag_Zt-19.0)^1.0),"starssize"]=3
}
#########################
# Assign object classes #
#########################
datafile0[,"class"]="ambiguous"
datafile0[starscol<0.5,"class"]="galaxy"
datafile0[starscol>2.5,"class"]="star"
datafile0[starssize<0.5,"class"]="galaxy"
datafile0[starssize>2.5,"class"]="star"
#########################################################################################
# Identify artifacts (no detection outside r+Z, too small and bright improbably colour) #
#########################################################################################
medianskyRMS=quantile(datafile0[,skyRMS_mean],0.5)
lowsky=datafile0[skyRMS_mean < medianskyRMS,skyRMS_mean]
crudcut=medianskyRMS+5*(medianskyRMS-quantile(lowsky,0.33))
#
datafile0[skyRMS_mean > crudcut,"class"]="artefact"
if (is.na(datafile0$mag_ut[1])==TRUE & is.na(datafile0$mag_gt[1])==TRUE & is.na(datafile0$mag_rt[1])==TRUE & is.na(datafile0$mag_it[1])==TRUE) {
	datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kt) > 2 & noIR < 1,"class"]="artefact"
	datafile0[log10(R50) < -0.4,"class"]="artefact"
	datafile0[mag_Zt-mag_Yt < -1.5 | mag_Zt-mag_Yt > 1.5,"class"]="artifact"
} else {
	datafile0[is.na(mag_gt)+is.na(mag_rt)+is.na(mag_it) > 2 & noOPT < 1,"class"]="artefact"
	datafile0[is.na(mag_Zt)+is.na(mag_Yt)+is.na(mag_Jt)+is.na(mag_Ht)+is.na(mag_Kt) > 2 & noIR < 1,"class"]="artefact"
	datafile0[log10(R50) < -0.4,"class"]="artefact"
	datafile0[(mag_rt-mag_Zt) < -0.75 & noIR < 1 & noOPT < 1,"class"]="artefact"
}
###############################
# Add in prior info from GAMA #
###############################
#radius=20.0-gamaprior$MAG
#radius=ifelse(radius<2.5,2.5,radius)
#gamamatch=coordmatch(coordcompare=datafile0[,list(RAmax,Decmax)],coordref=gamaprior[,list(RA,DEC)],rad=radius,inunitref="deg",inunitcompare="deg",radunit="asec")
#datafile0[gamamatch$bestmatch$compareID,"CATAID"]=gamaprior[gamamatch$bestmatch$refID,"CATAID"]
#datafile0[gamamatch$bestmatch$compareID,"Z"]=gamaprior[gamamatch$bestmatch$refID,"Z"]
#datafile0[gamamatch$bestmatch$compareID,"NQ"]=gamaprior[gamamatch$bestmatch$refID,"NQ"]
#datafile0[gamamatch$bestmatch$compareID,"PROB"]=gamaprior[gamamatch$bestmatch$refID,"PROB"]
###############################################
# Override classification if quality redshift #
###############################################
datafile0$class=ifelse(datafile0$Z>0.002 & datafile0$NQ>2,"galaxy",datafile0$class)
datafile0$class=ifelse(datafile0$Z>-0.002 & datafile0$Z<0.002 & datafile0$NQ>2,"star",datafile0$class)
######################
# Identify star mask #
######################
gaia0=subset(gaia0,gaia0$phot_g_mean_mag<18)
radius=10^(1.6-0.15*gaia0$phot_g_mean_mag)
radius=ifelse(gaia0$phot_g_mean_mag<6.0,5.0119,radius)
#
maskedregion=coordmatch(coordcompare=datafile0[,list(RAGAIA,DecGAIA)],coordref=gaia0[,list(ra,dec)],rad=radius,inunitref="deg",inunitcompare="deg",radunit="amin")
length(maskedregion$bestmatch$refID)
maskedobjs=as.numeric(maskedregion$ID[,])
maskedobjs=maskedobjs[maskedobjs>0]
length(maskedobjs)
datafile0[maskedobjs,"starmask"]=1
#############################
# select fragmented systems #
#############################
groups_merge = everything$pro_detect$group$groupsegID[everything$pro_detect$group$groupsegID$Ngroup>5,"groupID"]
#######################
# build new catalogue #
#######################
new_tot =  profoundCatMerge(everything$cat_tot, everything$cat_grp, everything$pro_detect$group$groupsegID, groups_merge)
new_col = everything$cat_col
new_seg = everything$pro_detect$segstats
new_col = new_col[match(new_tot$segID,new_col$segID),]
new_seg = new_seg[match(new_tot$segID,new_seg$segID),]
fdatafile <- as.data.table(cbind(new_seg, new_col, new_tot))

if (is.na(datafile0$mag_ut[1])==TRUE & is.na(datafile0$mag_gt[1])==TRUE & is.na(datafile0$mag_rt[1])==TRUE & is.na(datafile0$mag_it[1])==TRUE) {
	fdatafile0=as.data.table(fdatafile[,list(uniqueID,segID,xmax,ymax,xcen,ycen,RAcen,Deccen,RAmax,Decmax,flux_Zt,flux_err_Zt,flux_Yt,flux_err_Yt,flux_Jt,flux_err_Jt,flux_Ht,flux_err_Ht,flux_Kt,flux_err_Kt,flux_Zc,flux_err_Zc,flux_Yc,flux_err_Yc,flux_Jc,flux_err_Jc,flux_Hc,flux_err_Hc,flux_Kc,flux_err_Kc,sky_mean,skyRMS_mean,SB_N50_Zt,SB_N100_Zt,R50,R90,R100,N100,mag_Zc,mag_Zt,mag_Yc,mag_Yt,mag_Jc,mag_Jt,mag_Hc,mag_Ht,mag_Kc,mag_Kt)])
	} else {
fdatafile0=as.data.table(fdatafile[,list(uniqueID,segID,xmax,ymax,xcen,ycen,RAcen,Deccen,RAmax,Decmax,flux_ut,flux_err_ut,flux_gt,flux_err_gt,flux_rt,flux_err_rt,flux_it,flux_err_it,flux_Zt,flux_err_Zt,flux_Yt,flux_err_Yt,flux_Jt,flux_err_Jt,flux_Ht,flux_err_Ht,flux_Kt,flux_err_Kt,flux_uc,flux_err_uc,flux_gc,flux_err_gc,flux_rc,flux_err_rc,flux_ic,flux_err_ic,flux_Zc,flux_err_Zc,flux_Yc,flux_err_Yc,flux_Jc,flux_err_Jc,flux_Hc,flux_err_Hc,flux_Kc,flux_err_Kc,sky_mean,skyRMS_mean,SB_N50_Zt,SB_N100_Zt,SB_N50_rt,SB_N100_rt,R50,R90,R100,N100,mag_uc,mag_ut,mag_gc,mag_gt,mag_rc,mag_rt,mag_ic,mag_it,mag_Zc,mag_Zt,mag_Yc,mag_Yt,mag_Jc,mag_Jt,mag_Hc,mag_Ht,mag_Kc,mag_Kt)])
}
#
iauname=IAUID(fdatafile0$RAcen,fdatafile0$Deccen,name = 'WAVES',epoch = '_N_')
frameid=ceiling(min(fdatafile0$RAcen))*100.0+ceiling(min(fdatafile0$Deccen))
uberid=frameid*1e10+fdatafile0$uniqueID
seeing=log10(median(fdatafile0[mag_Zt > 16.0 & mag_Zt < 17 & R50 < 2.0,R50]))
fdatafile0=cbind(iauname,frameid,uberid,fdatafile0,seeing)
names(fdatafile0)[1] <- "IAUID"
names(fdatafile0)[2] <- "FrameID"
names(fdatafile0)[3] <- "uberID"
names(fdatafile0)[length(fdatafile0)] <- "seeing"
################################################################
# Determine extinction (E(B-V) values using Planck E(B-V) map) #
################################################################
waves=coordmatch(coordref=fdatafile0[,list(RAmax,Decmax)],coordcompare=ebv0[,list(RA_J2000,DEC_J2000)],rad=150,inunitref="deg",inunitcompare="deg",radunit="asec")
fdatafile0=cbind(fdatafile0,ebv0[waves$ID[,1],"EBV"])
#
fdatafile0 = as.data.table(cbind(fdatafile0,censep=as.numeric(sqrt((11000/2.0-fdatafile0[,"xmax"])^2+(11000/2.0-fdatafile0[,"ymax"])^2)),RAGAIA=fdatafile0$RAcen-gaiaraoff,DecGAIA=fdatafile0$Deccen-gaiadecoff,class="notclassified",duplicate=1,starscol=0,starssize=0,mask=0,noOPT=0,noIR=0,regroup=0,starmask=0,CATAID=0,Z=-999,NQ=0,PROB=0))
#######################################################################################
# Re-Correct all mags for attenuation using coefficients from Galactic Extinction DMU #
#######################################################################################
optmult=10^{0.4*(8.9-0.0)}
nirmult=10^{0.4*(8.9-30.0)}
#  
Au=5.155
Ag=3.793
Ar=2.751
Ai=2.086
AZ=1.71
AY=1.31
AJ=0.928
AH=0.592
AK=0.388
#
fdatafile0$mag_ut=ifelse(fdatafile0$mag_ut==-999,-999,fdatafile0$mag_ut-Au*fdatafile0$EBV)
fdatafile0$mag_gt=ifelse(fdatafile0$mag_gt==-999,-999,fdatafile0$mag_gt-Ag*fdatafile0$EBV)
fdatafile0$mag_rt=ifelse(fdatafile0$mag_rt==-999,-999,fdatafile0$mag_rt-Ar*fdatafile0$EBV)
fdatafile0$mag_it=ifelse(fdatafile0$mag_it==-999,-999,fdatafile0$mag_it-Ai*fdatafile0$EBV)
fdatafile0$mag_Zt=ifelse(fdatafile0$mag_Zt==-999,-999,fdatafile0$mag_Zt-AZ*fdatafile0$EBV)
fdatafile0$mag_Yt=ifelse(fdatafile0$mag_Yt==-999,-999,fdatafile0$mag_Yt-AY*fdatafile0$EBV)
fdatafile0$mag_Jt=ifelse(fdatafile0$mag_Jt==-999,-999,fdatafile0$mag_Jt-AJ*fdatafile0$EBV)
fdatafile0$mag_Ht=ifelse(fdatafile0$mag_Ht==-999,-999,fdatafile0$mag_Ht-AH*fdatafile0$EBV)
fdatafile0$mag_Kt=ifelse(fdatafile0$mag_Kt==-999,-999,fdatafile0$mag_Kt-AK*fdatafile0$EBV)
#
fdatafile0$mag_uc=ifelse(fdatafile0$mag_uc==-999,-999,fdatafile0$mag_uc-Au*fdatafile0$EBV)
fdatafile0$mag_gc=ifelse(fdatafile0$mag_gc==-999,-999,fdatafile0$mag_gc-Ag*fdatafile0$EBV)
fdatafile0$mag_rc=ifelse(fdatafile0$mag_rc==-999,-999,fdatafile0$mag_rc-Ar*fdatafile0$EBV)
fdatafile0$mag_ic=ifelse(fdatafile0$mag_ic==-999,-999,fdatafile0$mag_ic-Ai*fdatafile0$EBV)
fdatafile0$mag_Zc=ifelse(fdatafile0$mag_Zc==-999,-999,fdatafile0$mag_Zc-AZ*fdatafile0$EBV)
fdatafile0$mag_Yc=ifelse(fdatafile0$mag_Yc==-999,-999,fdatafile0$mag_Yc-AY*fdatafile0$EBV)
fdatafile0$mag_Jc=ifelse(fdatafile0$mag_Jc==-999,-999,fdatafile0$mag_Jc-AJ*fdatafile0$EBV)
fdatafile0$mag_Hc=ifelse(fdatafile0$mag_Hc==-999,-999,fdatafile0$mag_Hc-AH*fdatafile0$EBV)
fdatafile0$mag_Kc=ifelse(fdatafile0$mag_Kc==-999,-999,fdatafile0$mag_Kc-AK*fdatafile0$EBV)
#
fdatafile0$flux_ut=ifelse(fdatafile0$flux_ut==0.0,-999,fdatafile0$flux_ut*10^(0.4*(Au*fdatafile0$EBV))*optmult)
fdatafile0$flux_gt=ifelse(fdatafile0$flux_gt==0.0,-999,fdatafile0$flux_gt*10^(0.4*(Ag*fdatafile0$EBV))*optmult)
fdatafile0$flux_rt=ifelse(fdatafile0$flux_rt==0.0,-999,fdatafile0$flux_rt*10^(0.4*(Ar*fdatafile0$EBV))*optmult)
fdatafile0$flux_it=ifelse(fdatafile0$flux_it==0.0,-999,fdatafile0$flux_it*10^(0.4*(Ai*fdatafile0$EBV))*optmult)
fdatafile0$flux_Zt=ifelse(fdatafile0$flux_Zt==0.0,-999,fdatafile0$flux_Zt*10^(0.4*(AZ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Yt=ifelse(fdatafile0$flux_Yt==0.0,-999,fdatafile0$flux_Yt*10^(0.4*(AY*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Jt=ifelse(fdatafile0$flux_Jt==0.0,-999,fdatafile0$flux_Jt*10^(0.4*(AJ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Ht=ifelse(fdatafile0$flux_Ht==0.0,-999,fdatafile0$flux_Ht*10^(0.4*(AH*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Kt=ifelse(fdatafile0$flux_Kt==0.0,-999,fdatafile0$flux_Kt*10^(0.4*(AK*fdatafile0$EBV))*nirmult)
#
fdatafile0$flux_err_ut=ifelse(fdatafile0$flux_err_ut==0.0,-999,fdatafile0$flux_err_ut*10^(0.4*(Au*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_gt=ifelse(fdatafile0$flux_err_gt==0.0,-999,fdatafile0$flux_err_gt*10^(0.4*(Ag*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_rt=ifelse(fdatafile0$flux_err_rt==0.0,-999,fdatafile0$flux_err_rt*10^(0.4*(Ar*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_it=ifelse(fdatafile0$flux_err_it==0.0,-999,fdatafile0$flux_err_it*10^(0.4*(Ai*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_Zt=ifelse(fdatafile0$flux_err_Zt==0.0,-999,fdatafile0$flux_err_Zt*10^(0.4*(AZ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Yt=ifelse(fdatafile0$flux_err_Yt==0.0,-999,fdatafile0$flux_err_Yt*10^(0.4*(AY*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Jt=ifelse(fdatafile0$flux_err_Jt==0.0,-999,fdatafile0$flux_err_Jt*10^(0.4*(AJ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Ht=ifelse(fdatafile0$flux_err_Ht==0.0,-999,fdatafile0$flux_err_Ht*10^(0.4*(AH*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Kt=ifelse(fdatafile0$flux_err_Kt==0.0,-999,fdatafile0$flux_err_Kt*10^(0.4*(AK*fdatafile0$EBV))*nirmult)
#
fdatafile0$flux_uc=ifelse(fdatafile0$flux_uc==0.0,-999,fdatafile0$flux_uc*10^(0.4*(Au*fdatafile0$EBV))*optmult)
fdatafile0$flux_gc=ifelse(fdatafile0$flux_gc==0.0,-999,fdatafile0$flux_gc*10^(0.4*(Ag*fdatafile0$EBV))*optmult)
fdatafile0$flux_rc=ifelse(fdatafile0$flux_rc==0.0,-999,fdatafile0$flux_rc*10^(0.4*(Ar*fdatafile0$EBV))*optmult)
fdatafile0$flux_ic=ifelse(fdatafile0$flux_ic==0.0,-999,fdatafile0$flux_ic*10^(0.4*(Ai*fdatafile0$EBV))*optmult)
fdatafile0$flux_Zc=ifelse(fdatafile0$flux_Zc==0.0,-999,fdatafile0$flux_Zc*10^(0.4*(AZ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Yc=ifelse(fdatafile0$flux_Yc==0.0,-999,fdatafile0$flux_Yc*10^(0.4*(AY*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Jc=ifelse(fdatafile0$flux_Jc==0.0,-999,fdatafile0$flux_Jc*10^(0.4*(AJ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Hc=ifelse(fdatafile0$flux_Hc==0.0,-999,fdatafile0$flux_Hc*10^(0.4*(AH*fdatafile0$EBV))*nirmult)
fdatafile0$flux_Kc=ifelse(fdatafile0$flux_Kc==0.0,-999,fdatafile0$flux_Kc*10^(0.4*(AK*fdatafile0$EBV))*nirmult)
#
fdatafile0$flux_err_uc=ifelse(fdatafile0$flux_err_uc==0.0,-999,fdatafile0$flux_err_uc*10^(0.4*(Au*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_gc=ifelse(fdatafile0$flux_err_gc==0.0,-999,fdatafile0$flux_err_gc*10^(0.4*(Ag*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_rc=ifelse(fdatafile0$flux_err_rc==0.0,-999,fdatafile0$flux_err_rc*10^(0.4*(Ar*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_ic=ifelse(fdatafile0$flux_err_ic==0.0,-999,fdatafile0$flux_err_ic*10^(0.4*(Ai*fdatafile0$EBV))*optmult)
fdatafile0$flux_err_Zc=ifelse(fdatafile0$flux_err_Zc==0.0,-999,fdatafile0$flux_err_Zc*10^(0.4*(AZ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Yc=ifelse(fdatafile0$flux_err_Yc==0.0,-999,fdatafile0$flux_err_Yc*10^(0.4*(AY*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Jc=ifelse(fdatafile0$flux_err_Jc==0.0,-999,fdatafile0$flux_err_Jc*10^(0.4*(AJ*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Hc=ifelse(fdatafile0$flux_err_Hc==0.0,-999,fdatafile0$flux_err_Hc*10^(0.4*(AH*fdatafile0$EBV))*nirmult)
fdatafile0$flux_err_Kc=ifelse(fdatafile0$flux_err_Kc==0.0,-999,fdatafile0$flux_err_Kc*10^(0.4*(AK*fdatafile0$EBV))*nirmult)
#############################
# copy over derived columns #
#############################
fdatafile0$EBV = datafile0[match(fdatafile0$segID,datafile0$segID),EBV]
fdatafile0$class = datafile0[match(fdatafile0$segID,datafile0$segID),class]
fdatafile0$duplicate = datafile0[match(fdatafile0$segID,datafile0$segID),duplicate]
fdatafile0$starscol = datafile0[match(fdatafile0$segID,datafile0$segID),starscol]
fdatafile0$starssize = datafile0[match(fdatafile0$segID,datafile0$segID),starssize]
fdatafile0$mask = datafile0[match(fdatafile0$segID,datafile0$segID),mask]
fdatafile0$noOPT = datafile0[match(fdatafile0$segID,datafile0$segID),noOPT]
fdatafile0$noIR = datafile0[match(fdatafile0$segID,datafile0$segID),noIR]
fdatafile0$starmask = datafile0[match(fdatafile0$segID,datafile0$segID),starmask]
fdatafile0$CATAID = datafile0[match(fdatafile0$segID,datafile0$segID),CATAID]
fdatafile0$Z = datafile0[match(fdatafile0$segID,datafile0$segID),Z]
fdatafile0$NQ = datafile0[match(fdatafile0$segID,datafile0$segID),NQ]
fdatafile0$PROB = datafile0[match(fdatafile0$segID,datafile0$segID),PROB]
#################################
# set regroup flag if regrouped #
#################################
fdatafile0$regroup[fdatafile0$segID %in% groups_merge]=1
###################
# Omit NA columns #
###################
fdatafile0 = fdatafile0 %>% select_if(~sum(!is.na(.)) > 0)
#######################################
# select key columns to carry forward #
#######################################
trim=as.list("trim")
trim$rawcat <- datafile0
trim$finalcat <- fdatafile0
trim$groupcat <- everything$pro_detect$group$groupsegID
trim$detect_segim=everything$pro_detect$segim_orig
trim$dilate_segim=everything$pro_detect$segim
trim$final_segim=profoundSegimKeep(everything$pro_detect$segim,everything$pro_detect$group$groupim,groups_merge)
trim$header=everything$pro_detect$header
#########################
# save trimmed RDS file #
#########################
out = '/group/pawsey0160/waves/trim/'
saveRDS(trim,file=paste0(out,"RA",i,"DEC",j,".rds",sep=''))
#############
# Write CSV #
#############
out_csv = '/group/pawsey0160/waves/trim/'
write.csv(trim$finalcat,file=paste0(out_csv,"RA",i,"DEC",j,".csv",sep=''),na="-999",row.names=FALSE)
#


		} # for RA
	} # for DEC
#############
# Make plot #
#############
png(paste0(stub,"swarpit/plots",region,"png"),width=25.0,height=25.0,units="cm",res=240)
#
magimageWCS(everything$pro_detect$image, header=everything$pro_detect$header)
points(fdatafile0[starmask > 0,list(xmax,ymax)],pch=16,col="orange",cex=0.25)
points(fdatafile0[CATAID > 0 & Z > 0.002 & NQ > 2 & starmask < 1,list(xmax,ymax)],pch=".",col="green",cex=0.25)
points(fdatafile0[class=="artefact",list(xmax,ymax)],pch=4,col="cyan",cex=0.25)
points(fdatafile0[noIR > 0,list(xmax,ymax)],pch=18,col="magenta",cex=0.5)
points(fdatafile0[noOPT >0,list(xmax,ymax)],pch=18,col="yellow",cex=0.5)
points(fdatafile0[CATAID < 1 & mag_rt < 19.8 & mag_rt > -999 & starmask < 1 & class=="galaxy",list(xmax,ymax)],pch=5,col="red",cex=0.75)
title(main=region,sub=paste0("GAMA galaxies:",length(fdatafile0[CATAID >0 & Z > 0.002 & NQ > 2,segID]),"     Artefacts: ",length(fdatafile0[class=="artefact",segID]),"     Not in GAMA: ",length(fdatafile0[CATAID < 1 & mag_rt < 19.8 & mag_rt > 0.0 & starmask < 1 & class=="galaxy",segID]),"     Masked: ",length(fdatafile0[starmask > 0,segID]),"     NoOPT: ",length(fdatafile0[noOPT > 0,segID]),"     NoIR: ",length(fdatafile0[noIR > 0,segID])))
#
dev.off()
