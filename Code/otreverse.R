otreverse = function(ldiff){
  
  ## Inputs
  # ldiff: Px1 matrix of forecast loss differences from two models
  
  # Length of forecast evaluation window
  P = nrow(ldiff)
  
  # Lower and upper trimming values for test statistic
  lower = round(0.15*P)
  upper = round(0.85*P)
  
  # Bandwidth of HAC estimator
  bw = P^(1/3)
  
  # HAC Variance of forecast loss differential
  ldiffvar = var(ldiff[,1])
  #   ldifftilde = ldiff-mean(ldiff)
  #   for (j in 1:bw){
  #     gammatemp = t(ldifftilde[(j+1):P,1])%*%ldifftilde[1:(P-j),1]/P
  #     ldiffvar = ldiffvar+2*(1-j/(bw+1))*gammatemp
  #   }
  
  # One-time reversal test statistic statistic
  lm1 = sum(ldiff[,1])^2/(ldiffvar*P)
  lm2 = matrix(0,upper-lower+1,1)
  j = 1
  for (t in lower:upper){
    lm2[j,1] = 1/(ldiffvar*t*(1-t/P))*(sum(ldiff[1:t,1])-(t/P)*sum(ldiff[,1]))^2
    j = j+1
  }
  otrevstat = max(lm1+lm2)
  
  # Test decision
  rej = 0
  if (otrevstat>9.83){
    rej = 1
  }
  
  # Causes of rejection
  rejlm1 = 0 # One model is constantly better than the other
  rejlm2 = 0 # There is a one-time reversal in relative forecast performance
  if (rej==1){
    if (lm1>3.84){
      rejlm1 = 1
    }
    if (max(lm2)>8.85){
      rejlm2 = 1
    }
  }
  
  # Estimation of break date
  breakdate = NA
  plotperform = matrix(NA,P)
  if (rejlm2==1){
    breakdate = which.max(lm2)+lower-1
    plotperform[1:breakdate] = sum(ldiff[1:breakdate,1])/breakdate
    plotperform[(breakdate+1):P] = sum(ldiff[(breakdate+1):P,1])/(P-breakdate)
  }
  
  results = list(otrevstat=otrevstat,rej=rej,rejlm1=rejlm1,rejlm2=rejlm2,breakdate=breakdate,plotperform=plotperform)
  return(results)
}