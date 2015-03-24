WS <- c(ls())
rm(list=WS)
load('h:/git/bubblesbreakdowns/results/IPT hor 13 rolling_ARX.RData')
load('h:/git/bubblesbreakdowns/results/IPT hor 13 rolling_AR.RData')

## Forecast evaluation
yeval = set.rt$IPT
yeval = as.matrix(1200*log(yeval[(1+horizon):length(yeval)]/yeval[1:(length(yeval)-horizon)])/horizon,nrow=1)
yeval = yeval[!is.na(yeval)]

# Remove last #horizon forecasts for which no evaluation is possible
vint.dates = vint.dates[1:(length(vint.dates)-horizon-1)]
Neval = length(vint.dates)

forecast.ar = as.numeric(forecast.ar[1:Neval,,drop=F])
forecast.arx = forecast.arx[,1:Neval,drop=F]
msr.ar = as.numeric(msr.ar[1:Neval,,drop=F])
msr.arx = msr.arx[,1:Neval,drop=F]

forecast.arx = rbind(t(forecast.ar),forecast.arx)
row.names(forecast.arx)[1]="AR"
msr.arx = rbind(t(msr.ar),msr.arx)
row.names(msr.arx)[1]="AR"
Nmodels = nrow(forecast.arx)

yeval = as.matrix(yeval[(length(yeval)-Neval+1):length(yeval)],nrow=1)

# Timeline of forecast margins
timeline=seq((1983+7/12),(2013+9/12),by=(1/12))

# Forecast errors
error.arx = forecast.arx-matrix(1,Nmodels,1)%*%t(yeval)
sqerror.arx = error.arx^2

# Surprise loss
SL.arx = sqerror.arx-msr.arx

# Variances
sqerror.arx.demeaned=sqerror.arx-rowMeans(sqerror.arx)%*%matrix(1,nrow=1,ncol=Neval)
SLL.arx=cov(t(sqerror.arx.demeaned))
SLL.arx=diag(SLL.arx)

# Parameters
lambda=2/3*(max.obs/Neval)

# HAC variance estimator of demeaned surprise losses
bw=0#floor(Neval^(1/3))# rounded down
SLL.arx = matrix(0,nrow=1,ncol=Nmodels)
for (n in 1:Nmodels){
  hacest<-function(sqerror.demeaned,Neval,bw){
    SLL = 0
    for (j in 1:bw){
      Gamma = t(sqerror.demeaned[(1+j):Neval])%*%sqerror.demeaned[1:(Neval-j)]/Neval
      SLL = SLL+2*(1-j/(bw+1))*Gamma
    }
    SLL = SLL+t(sqerror.demeaned)%*%sqerror.demeaned/Neval
  }
  SLL.arx=apply(sqerror.arx.demeaned,1,hacest,Neval,bw)
}
# Variance estimator out-of-sample losses
sigma2 = lambda*SLL.arx

# regression of SL on themselves
pmax = 12
source(paste(DirCode,'/olsself.R',sep=''))
surprise.BIC = matrix(NA,nrow=pmax,ncol=Nmodels)
for (p in 1:pmax){
  surprise.res=apply(SL.arx,1,olsself,p)
  surprise.BIC[p,]=sapply(surprise.res,function(x) x$BIC)
}
surprise.pstar = apply(surprise.BIC, 2, which.min)
for (n in 1:Nmodels){
  surprise.res[[n]]=olsself(SL.arx[n,],surprise.pstar[n])
}

surprise.coef=sapply(surprise.res,function(x) x$b)
surprise.regressors=sapply(surprise.res,function(x) x$Z)
surprise.resid=sapply(surprise.res,function(x) x$res)
surprise.fit=sapply(surprise.res,function(x) x$yfit)

# Confidence Interval
source(paste(DirCode,'/FB_Wald_CI.R',sep=''))
surprise.CI = list()
for (n in 1:Nmodels){
  trash = FB_Wald_CI(surprise.res[[n]]$Z,surprise.res[[n]]$res,surprise.res[[n]]$b,Neval,max.obs,sigma2[n],1,"rolling",bw,0.05)
  surprise.CI[[n]]=trash$SLfitCI
}

## Plots
dataplot = cbind(matrix(0,nrow=Neval-surprise.pstar[1],ncol=1),surprise.fit[[1]],surprise.CI[[1]])
matplot(timeline[(surprise.pstar[1]+1):Neval],dataplot,type="l")

## Extract data for plots
SL.arx.fit = matrix(NA,nrow=Nmodels,ncol=Neval)
SL.arx.CI = matrix(NA,nrow=Nmodels,ncol=Neval)
for (n in 1:Nmodels){
  trashfit = surprise.fit[[n]]
  trashCI = surprise.CI[[n]]
  SL.arx.fit[n,] = rbind(matrix(NA,nrow=surprise.pstar[n],ncol=1),trashfit)
  SL.arx.CI[n,] = rbind(matrix(NA,nrow=surprise.pstar[n],ncol=1),trashCI)
}

# Best model by rank
fcastranks = apply(sqerror.arx,2,rank)
fcastranksavg = apply(fcastranks,1,mean)
# Best model by rank
which.min(fcastranksavg)
# Forecast breakdowns
modbreak=(SL.arx.CI>0)*1
Nmodbreak=rowSums(modbreak,na.rm=T)
cor.test(fcastranksavg,Nmodbreak)

write.table(SL.arx,file="h:/git/bubblesbreakdowns/results/SLarx.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
write.table(SL.arx.fit,file="h:/git/bubblesbreakdowns/results/SLarxfit.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
write.table(SL.arx.CI,file="h:/git/bubblesbreakdowns/results/SLarxCI.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
write.table(modbreak,file="h:/git/bubblesbreakdowns/results/breakdowns.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
write.table(timeline,file="h:/git/bubblesbreakdowns/results/timeline.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)


# # regression of SL on themselves and on bubble indicators
# bubble=read.table('H:/git/bubblesbreakdowns/Results/stockhousebubble.csv')
# pmax = 12
# qmax = 1
# source(paste(DirCode,'/olsarx.R',sep=''))
# surprise.pstar = matrix(NA,nrow=1,ncol=Nmodels)
# for (n in 1:Nmodels){
#   BIC = matrix(NA,nrow=pmax,ncol=qmax)
#   for (p in 1:pmax){
#     trash = olsarx(SL.arx[n,],bubble,p,qmax)
#     BIC[p,] = trash$BIC
#     surprise.pstar[n] = which.min(BIC)
#   }
# }
# for (n in 1:Nmodels){
#   surprise.res[[n]]=olsarx(SL.arx[n,],bubble,surprise.pstar[n],qmax)
# }
# 
# surprise.coef=sapply(surprise.res,function(x) x$b)
# surprise.regressors=sapply(surprise.res,function(x) x$Z)
# surprise.resid=sapply(surprise.res,function(x) x$res)
# surprise.fit=sapply(surprise.res,function(x) x$yfit)
# 
# # Confidence Interval
# source(paste(DirCode,'/FB_Wald_CI.R',sep=''))
# surprise.CI = list()
# for (n in 1:Nmodels){
#   trash = FB_Wald_CI(surprise.res[[n]]$Z,surprise.res[[n]]$res,surprise.res[[n]]$b,Neval,max.obs,sigma2[n],1,"rolling",0,0.05)
#   surprise.CI[[n]]=trash$SLfitCI
# }
# 
# ## Plots
# dataplot = cbind(matrix(0,nrow=Neval-surprise.pstar[1],ncol=1),surprise.fit[[1]],surprise.CI[[1]])
# matplot(timeline[(surprise.pstar[1]+1):Neval],dataplot,type="l")
# 
# ## Extract data for plots
# SL.arx.fit = matrix(NA,nrow=Nmodels,ncol=Neval)
# SL.arx.CI = matrix(NA,nrow=Nmodels,ncol=Neval)
# 
# for (n in 1:Nmodels){
#   trashfit = surprise.fit[[n]]
#   trashCI = surprise.CI[[n]]
#   SL.arx.fit[n,] = rbind(matrix(NA,nrow=surprise.pstar[n],ncol=1),trashfit)
#   SL.arx.CI[n,] = rbind(matrix(NA,nrow=surprise.pstar[n],ncol=1),trashCI)
# }
# 
# write.table(SL.arx,file="h:/git/bubblesbreakdowns/results/SLarxSH.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
# write.table(SL.arx.fit,file="h:/git/bubblesbreakdowns/results/SLarxfitSH.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
# write.table(SL.arx.CI,file="h:/git/bubblesbreakdowns/results/SLarxCISH.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)
# write.table(timeline,file="h:/git/bubblesbreakdowns/results/timelineSH.csv",quote=T,append=F,sep=",",eol = "\n", na = "NaN", dec = ".", row.names = F,col.names = F)