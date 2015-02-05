# target variable
horizon.nonr=12
df.unrevised.nobs=nrow(df.unrevised)
y=df.unrevised[,'PCPIX',drop=F]
y=chg(y,12)
vint.dates.obs=which(row.names(y)%in%vint.dates)
vint.dates.obs=vint.dates.obs+horizon.nonr
vint.dates.obs=vint.dates.obs[vint.dates.obs<df.unrevised.nobs]
y=y[vint.dates.obs,1,drop=F]
# forecast values

forecast.value=sapply(forecast.all[1:184],function(x) x[,'bma.fc'] )
forecast.value=t(forecast.value)
row.names(forecast.value)=vint.dates[1:184]
model.names=forecast.all[[1]]$names
model.n=nrow(forecast.all[[1]])
colnames(forecast.value)=model.names

# forecast errors
y=y[1:184,1,drop=F]
forecast.error=forecast.value-as.matrix(y)%*%matrix(1,nrow=1,ncol=model.n)

# insample fits

forecast.residvar=t(sapply(forecast.all[1:184],function(x) x[,'residvar'] ))
row.names(forecast.residvar)=vint.dates[1:184]
colnames(forecast.residvar)=model.names

forecast.nobs=t(sapply(forecast.all[1:184],function(x) x[,'nobs'] ))
row.names(forecast.nobs)=vint.dates[1:184]
colnames(forecast.nobs)=model.names

forecast.sqerror=forecast.error^2
# surprise losses
surpriseloss=forecast.sqerror-forecast.residvar

# 
forecast.sqerror.demeaned=forecast.sqerror-matrix(1,nrow=184,ncol=1)%*%colMeans(forecast.sqerror)
SLL=cov(forecast.sqerror.demeaned)
SLL=diag(SLL)

# for lambda, we might need to fix the number of insample observations
min(apply(forecast.nobs,1,min)) #minimum nobs
# if nobs per model(!) is constant - more or less - it is ok.
surprise.m=round(colMeans(forecast.nobs),0)

# here, only for rolling regression.
surprise.lambda=2/3*(surprise.m/184)

# lambda=1
# if (s=='fixed'){lambda = 1+n/m}
# if (s=='rolling'& n<m){lambda = 1-(1/3)*(n/m)^2}
# if (s=='rolling'&n>=m){lambda = (2/3)*(m/n)}


# HAC variance estimator of demeaned surprise losses
nvint=184
bw=ceiling(184^(1/3))# rounded up, to be save (???check)
hacest<-function(forecast.sqerror.demeaned,nvint,bw){
        SLL = 0
     for (j in 1:bw){
             Gamma = t(forecast.sqerror.demeaned[(1+j):nvint])%*%forecast.sqerror.demeaned[1:(nvint-j)]/nvint
             SLL = SLL+2*(1-j/(bw+1))*Gamma
     }
     SLL = SLL+t(forecast.sqerror.demeaned)%*%forecast.sqerror.demeaned/nvint
}
SLL=apply(forecast.sqerror.demeaned,2,hacest,nvint,bw)
# Variance estimator out-of-sample losses
sigma2 = surprise.lambda*SLL

# regression of SL on themselves
pstar=3
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


surprise.wald=vector('list',model.n)
for (i in 1:model.n){
        surprise.wald[[i]]=FB_Wald_CI(surprise.res[[i]]$Z,surprise.res[[i]]$res,surprise.res[[i]]$b,nvint,surprise.res[[i]]$m,
                surprise.res[[i]]$sigma2,scheme='rolling',bw=bw)
}


surprise.ci=sapply(surprise.wald,function(x) x$SLfitCI)



