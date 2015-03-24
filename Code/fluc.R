fluc = function(ldiff,mu){
  
  ## Inputs
  # ldiff: Px1 matrix of forecast loss differences from two models
  # mu: Length of centered rolling window used for FLUC Test relative to forecast evaluation window
  
  # Length of forecast evaluation window
  P = nrow(ldiff)
  
  # Length of centered rolling window used for FLUC Test
  m = mu*P
  
  # Bandwidth of HAC estimator
  bw = P^(1/3)
  
  # HAC Variance of forecast loss differential
  ldiffvar = var(ldiff[,1])
#   ldifftilde = ldiff-mean(ldiff)
#   for (j in 1:bw){
#     gammatemp = t(ldifftilde[(j+1):P,1])%*%ldifftilde[1:(P-j),1]/P
#     ldiffvar = ldiffvar+2*(1-j/(bw+1))*gammatemp
#   }
  
  # Fluc Test statistic
  flucstat = matrix(0,P-m)
  i = 1
  for (j in (m/2+1):(P-m/2)){
    ldifftemp = ldiff[(j-m/2):(j+m/2),1]
    flucstat[i,1] = sum(ldifftemp)/sqrt(ldiffvar*m)
    i = i+1
  }
  
  # Critical values from Table 1 of Giacomini & Rossi
  critvals = c(3.39, 3.17, 3.02, 2.89, 2.77, 2.63, 2.56, 2.43, 2.24)
  critvalind = round(mu*10)
  crit = critvals[critvalind]*matrix(1,P-m)
  
  rej = 0
  if (any(abs(flucstat)>crit)){
    rej=1
  }
  
  results = list(flucstat=flucstat,crit=crit,rej=rej)
  return(results)
}