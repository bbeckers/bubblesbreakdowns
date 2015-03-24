olsarx=function(y,x,p,q){
## Header
# This function runs an OLS estimation of an autoregressive model with K
#
# Inputs:
# y: Tx1 time series of level variable y
# p: Number of lags to be included of difference of y
#
# Output:
# b: (p+1)x1 vector of estimated parameters
# res: (T-p)x1 vector of residuals
# BIC: Bayesian Information Criterion for p lags.
# Z: (T-p)x(p+1) matrix of regressors

## Function
y=as.matrix(y)
T = nrow(y)

Zp = matrix(NA,T-p,p)
for (i in 1:p){
Zp[,i] = y[(p+1-i):(T-i)]
}
Tp = T-p
if (q>0){
  Zq = matrix(NA,T-q,q)
  for (i in 1:q){
    Zq[,i] = x[(q+1-i):(T-i),]
  }
  Tq = T-q
  if (p>q){
    Z = cbind(Zp,Zq[(p-q+1):Tq,])
  }
  if (p<q){
    Z = cbind(Zp[(q-p+1):Tp,],Zq)
  }
  if (p==q){
    Z = cbind(Zp,Zq)
  }
}else{
  Z = Zp
}
Tz = nrow(Z)

Z = cbind(matrix(1,Tz,1),Z)
y = y[(max(p,q)+1):T]

# OLS estimator
b = solve(t(Z)%*%Z)%*%t(Z)%*%y
yfit = Z%*%b
res = y-yfit
sigma_u = t(res)%*%res/(T-ncol(Z)-1)
# Information criteria for determining the optimal lag length
BIC = T*log(sigma_u)+(ncol(Z)+1)*log(T)

results = list(b=b,res=res,yfit=yfit,BIC=BIC,Z=Z)
return(results)
}