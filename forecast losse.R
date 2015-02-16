# forecast losse ----------------------------------------------------------
# This file first gets all the data used in the experiment, drops those models that are not useful,
# that is, that could not be computed for the whole horizon (zeros prevent change rate computation, insample
# fit is not appropriate at some stage), and collects the in- and out-of-sample error to construct the SL.
# Finally, the SL-test is conducted with the endogeneous as explanatory variable.
# ATTENTION: Target variable and horizon need to be set.
# setting some values

# Settings ----------------------------------------------------------------
horizon=12 # inflation is published with on month lag, but in the following, the variables
# are only lagged.
target='IPT' # core inflation
max.lag=12 # maximum lag length to be considered 
max.obs=-1 # maximum number of past observations to be considered (rolling estimation);

# Preparing raw data ------------------------------------------------------
# The result file to be loaded only contains the forecast.all list.
# Thus, first the raw data have to be loaded, as well. 
DirCode='C:/Users/Dirk/Documents/GitHub/bubblesbreakdowns'
DirCode='h:/Git/bubblesbreakdowns'

# sourcing necassary scripts for estimation and forecast (for a description see "olsbmalag.Rmd")
source(paste(DirCode,'/lag.exact.R',sep=''))
source(paste(DirCode,'/diff.R',sep=''))
source(paste(DirCode,'/chg.R',sep=''))
source(paste(DirCode,'/bmafo.R',sep=''))
library("BMA")

# setting some values
horizon=12 # inflation is published with on month lag, but in the following, the variables
# are only lagged.
target='IPT' # core inflation
max.lag=12 # maximum lag length to be considered 
max.obs=-1 # maximum number of past observations to be considered (rolling estimation);
# setting negative window turns that of (recursive estimation)
# rolling window size will be: max.obs-max.lag-horizon, as lags need to be considered for
# estimation. 

# loading realtime data sets and unrevised data 
load(paste(DirCode,'/data/sets.Rdata',sep=''))
df.unrevised=read.csv(paste(DirCode,'/data/unrevised data.csv'
                            ,sep='')
                      ,sep=','
                      ,na.strings='NaN'
                      ,row.names=1
)
overview.rt=read.csv(paste(DirCode,'/overview.csv',sep=''),row.names=1)
overview.nr=read.csv(paste(DirCode,'/data/NonrevData overview.csv',sep=''),sep=',')
overview.nr[,1]=gsub('-','.',overview.nr[,1])
# Saving a complete version of unrevised complete (necessary to get target variable without lags)
df.unrevised.compl=df.unrevised

# dropping the saisonally adjusted series for now
df.unrevised=df.unrevised[,-grep('.SA',colnames(df.unrevised))]
overview.nr=overview.nr[-grep('.SA',overview.nr[,1]),]

# dropping "MCOILWTICO" as it is too short (starting in 1986, first iteration here is 1983:7)
df.unrevised=df.unrevised[,-grep("MCOILWTICO",colnames(df.unrevised))]
overview.nr=overview.nr[-grep("MCOILWTICO",overview.nr[,1]),]

# transforming rownames (dates) of unrevised data accordingly
df.ur.rownames=row.names(df.unrevised)
dates=strsplit(df.ur.rownames,'\\.')
dates=sapply(dates,function(x) x[c(3,2)])
dates=apply(dates,2,function(x) paste(x,collapse=':'))
row.names(df.unrevised)=dates
df.unrevised=df.unrevised[,as.character(overview.nr[,1])]
variables.tlag=overview.nr[overview.nr$Publication.Lag!=0,'Abbreviation']
for (var.tlag in variables.tlag){
        df.unrevised[,var.tlag]=lag.exact(df.unrevised[,var.tlag,drop=F]
                                          ,overview.nr[overview.nr$Abbr==var.tlag
                                                       ,'Publication.Lag'])
}

# getting vintage dates that correspond to dates in rownames
vint.names=names(sets)
vint.last19=grep('99M12',vint.names)# last vintage of the 20th century
vint.dates=vint.names
vint.dates[1:vint.last19]=paste('19',vint.dates[1:vint.last19],sep='')
vint.dates[(vint.last19+1):length(vint.dates)]=paste('20',vint.dates[(vint.last19+1):length(vint.dates)],sep='')
vint.missing.zeros=grep('M[1-9]$',vint.dates)
# vint.dates[vint.missing.zeros]=paste('0',vint.dates[vint.missing.zeros],sep='')
vint.dates[vint.missing.zeros]=gsub('M',':0',vint.dates[vint.missing.zeros])
vint.dates=gsub('M',':',vint.dates)
vintage=data.frame(name=vint.names,date=vint.dates)

# dates in data matrices, both realtime and unrevised, must 
# at least contain vintages (observation 2015.1 not existing yet, but vintage 2015.1).
# If not a line needs to be added
aux.match=match(vintage$date,row.names(df.unrevised))
missing.dates=vintage$date[is.na(aux.match)]
missing.dates=as.character(missing.dates)
df.unrevised[missing.dates,]=NA

# are sets to short, too?
set.rt=sets[[193]]
aux.match=match(vintage$date,row.names(set.rt))
missing.dates=vintage$date[is.na(aux.match)]
missing.dates=as.character(missing.dates)
# add missing lines
sets=lapply(sets,function(x){x[missing.dates,]=NA
                             return(x)})


# Getting forecast all and eliminate useless models ------------------------------------------------------
# Some of the models have not been usefull at the forecast origin from an insample fit perspective.
# Furthermore, due to zeros, some change rate models have not been used over the whole number of vintages.

load(paste(DirCode,'/Results/',target,' hor ',horizon,' recursive.RData',sep=''))
# zero eliminated modelsK
N_mod_vint=sapply(forecast.all,function(x) nrow(x))
min.mod=which.min(N_mod_vint)
min.mod.names=as.character(forecast.all[[min.mod]]$names)
Nvint=length(forecast.all)
Notin=matrix(NA,Nvint,10)
for (vint in 1:Nvint){
        test=as.character(forecast.all[[vint]]$names)
        notin=!test%in%min.mod.names
        n=length(test[notin])
        if (n>0){
                Notin[vint,1:n]=test[notin]      
        }
}
# Notin contains all models that do not appear in all models
eliminate=as.vector(Notin)
eliminate=eliminate[is.na(eliminate)==F]
eliminate=unique(eliminate)

for (vint in 1:Nvint){
        del=forecast.all[[vint]]$names%in%eliminate
        forecast.all[[vint]]=forecast.all[[vint]][!del,]
        
}  
text=paste(eliminate,colapse='',sep='')

useful=sapply(forecast.all,function(x) x$useful)
useful.ind=rowSums(useful)==Nvint
useful.nvint=rowSums(useful)
use.ful.models=colSums(useful)
plot(use.ful.models,main='number of useful models\n each vintages',xlab='vintage',ylab='number of useful models',type='l')
hist(useful.nvint,main='histogram of how many\n vintages each model is useful',ylab='number of models',xlab='number of vintages when useful')

# dropping useless models after saving complete results.
forecast.all.s=forecast.all
# forecast.all=forecast.all.s
for (i in 1:length(forecast.all)){
        forecast.all[[i]]=forecast.all[[i]][useful.ind,]
}



# Starting the forecast analysis ------------------------------------------

# target variable: realized values
df.unrevised.nobs=nrow(df.unrevised)
y=df.unrevised[,'PCPIX',drop=F]
if (target=='IPT'){
        y=sets[[length(sets)]][,target,drop=F]
}

# na - treatment to compute change rates
y.isnan=which(is.na(y)==F)
y.isnan.short=y.isnan[(horizon+1):length(y.isnan)]
y=chg(y[y.isnan,,drop=F],horizon)

# observations which have not been forecasted need to be dropped
y=y[is.na(y)==F,,drop=F]
# the dates that correspond to the vintages
vint.dates.obs=which(row.names(y)%in%vint.dates)
# the dates that correspond to the forecasts
target.dates.obs=vint.dates.obs+horizon
# some of the dates are useless as there are still no forecasts available
if (target=='PCPIX'){
        target.dates.obs=target.dates.obs[target.dates.obs<df.unrevised.nobs]}
if (target=='IPT'){
        target.dates.obs=target.dates.obs[1:which(target.dates.obs==nrow(y))]
}
y=y[target.dates.obs,1,drop=F]
target.dates=row.names(y)

# number of the first N vintages that are useful as their realization is not in the future
Nvint.useful=nrow(y)
vint.dates=vint.dates[1:Nvint.useful]

# adjust forecast.all to the useful periods
for (i in (Nvint:(Nvint.useful+1))){forecast.all[[i]]=NULL}


forecast.value=sapply(forecast.all,function(x) x[,'bma.fc'] )
forecast.value=t(forecast.value)
row.names(forecast.value)=vint.dates
model.names=forecast.all[[1]][,'names']
model.n=length(model.names)
colnames(forecast.value)=model.names

# forecast errors
forecast.error=forecast.value-as.matrix(y)%*%matrix(1,nrow=1,ncol=model.n)

# insample fits
forecast.msr=t(sapply(forecast.all,function(x) x[,'msr'] ))
row.names(forecast.msr)=vint.dates
colnames(forecast.msr)=model.names

forecast.nobs=t(sapply(forecast.all,function(x) x[,'nobs'] ))
row.names(forecast.nobs)=vint.dates
colnames(forecast.nobs)=model.names

forecast.sqerror=forecast.error^2
# surprise losses
surpriseloss=forecast.sqerror-forecast.msr

# 
forecast.sqerror.demeaned=forecast.sqerror-matrix(1,nrow=Nvint.useful,ncol=1)%*%colMeans(forecast.sqerror)
SLL=cov(forecast.sqerror.demeaned)
SLL=diag(SLL)

# for lambda, we might need to fix the number of insample observations
min(apply(forecast.nobs,1,min)) #minimum nobs
# if nobs per model(!) is constant - more or less - it is ok.
surprise.m=round(colMeans(forecast.nobs),0)

# here, only for rolling regression.
surprise.lambda=2/3*(surprise.m/Nvint.useful)

# lambda=1
# if (s=='fixed'){lambda = 1+n/m}
# if (s=='rolling'& n<m){lambda = 1-(1/3)*(n/m)^2}
# if (s=='rolling'&n>=m){lambda = (2/3)*(m/n)}


# HAC variance estimator of demeaned surprise losses

bw=ceiling(Nvint.useful^(1/3))# rounded up, to be save (???check)
hacest<-function(forecast.sqerror.demeaned,Nvint.useful,bw){
        SLL = 0
     for (j in 1:bw){
             Gamma = t(forecast.sqerror.demeaned[(1+j):Nvint.useful])%*%forecast.sqerror.demeaned[1:(Nvint.useful-j)]/Nvint.useful
             SLL = SLL+2*(1-j/(bw+1))*Gamma
     }
     SLL = SLL+t(forecast.sqerror.demeaned)%*%forecast.sqerror.demeaned/Nvint.useful
}
SLL=apply(forecast.sqerror.demeaned,2,hacest,Nvint.useful,bw)
# Variance estimator out-of-sample losses
sigma2 = surprise.lambda*SLL

# regression of SL on themselves
pstar=3
source(paste(DirCode,'/olsself.R',sep=''))
surprise.res=apply(surpriseloss,2,olsself,pstar)

surprise.coef=sapply(surprise.res,function(x) x$b)
surprise.regressors=sapply(surprise.res,function(x) x$Z)
surprise.resid=sapply(surprise.res,function(x) x$res)

t1=surprise.regressors[,1]
t2=surprise.coef[,1]
t3=surprise.res[[1]]$Z %*% surprise.res[[1]]$b

surprise.fit=sapply(surprise.res,function(x) x$Z%*%x$b)
for (i in 1:model.n){
        surprise.res[[i]]$m=surprise.m[i]
        surprise.res[[i]]$sigma2=sigma2[i]
}

source(paste(DirCode,'/FB_Wald_CI.R',sep=''))
surprise.wald=vector('list',model.n)
for (i in 1:model.n){
        surprise.wald[[i]]=FB_Wald_CI(surprise.res[[i]]$Z,surprise.res[[i]]$res,surprise.res[[i]]$b,Nvint.useful,surprise.res[[i]]$m,
                surprise.res[[i]]$sigma2,scheme='rolling',bw=bw)
}


surprise.ci=sapply(surprise.wald,function(x) x$SLfitCI)



