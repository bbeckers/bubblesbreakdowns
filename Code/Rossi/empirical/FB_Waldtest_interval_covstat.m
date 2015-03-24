function interval = FB_Waldtest_interval_covstat(deltahat,Z,m,n,tau,scheme,hac_lag,R,reshat,condhomosk)

% Constructs confidence band for the Conditional Test of Giacomini and Rossi
% (2005) for the covariance stationary case
%
%
% OUTPUT: interval, an nx2 matrix containing fitted Surprise Losses in the
% first column, and the one-sided 90% confidence band in the second column
% INPUTS: deltahat, the OLS coefficient estimate in the regression of the surprise losses on explanatory
%         variables Z, written as a 1xr column vector
%         Z, the nxr matrix of explanatory variables (typically including a
%         column of 1)
%         m, the in-sample size
%         n, the out-of-sample size
%         tau, the forecast horizon
%         scheme, the forecasting scheme. scheme = 1 for fixed; scheme = 2
%                 for rolling; scheme = 3 for recursive
%         hac_lag, the truncation lag for the Newey-West estimator (a
%                  rule-of-thumb choice is n^(1/3)) of the asymptotic
%                  variance
%         R, the matrix of restrictions (dimension #restr*r)
%         reshat, the residuals of the regression
%         condhomosk=1 if conditional homoskedasticity else for conditional
%                  heteroskedasticity


T = m + n + tau - 1;
pi = n / m ; 
Znotdemeaned=Z;
Z=[Z(:,1),Z(:,2:cols(Z))- mmult(mean( Z(:,2:cols(Z)) ,1 )',ones(cols(Z)-1,rows(Z)))'  ];
WW=[1,- mean( Znotdemeaned(:,2:cols(Z)) ,1 )* inv((Z(:,2:cols(Z))'*Z(:,2:cols(Z)))/rows(Z)); zeros( cols(Z(:,2:cols(Z))),1  ) , inv((Z(:,2:cols(Z))'*Z(:,2:cols(Z)))/rows(Z))]; 

    % --------------------------------  FIXED SCHEME  --------------------
  if scheme == 1;
    lambda=1+n/m; lambdaLZL=1; 
    scaling=[lambda,ones(1,cols(Z)-1).*lambdaLZL;ones(cols(Z)-1,1).*lambdaLZL,ones(cols(Z)-1,cols(Z)-1)]; 
    samplevar=NeweyWest([reshat, mmult( reshat, Z(:,2:cols(Z)) ) ],hac_lag).*scaling;
    if condhomosk==1; samplevar=diag(diag(samplevar)); end; 
    omegahat =  WW*samplevar*WW';

    % --------------------------------  ROLLING SCHEME  --------------------
  elseif scheme == 2;
      if pi<1; 
        lambda=1-(1/3)*(n/m)^2;
        lambdaLZL=1-(n/m)/2; 
        scaling=[lambda,ones(1,cols(Z)-1).*lambdaLZL;ones(cols(Z)-1,1).*lambdaLZL,ones(cols(Z)-1,cols(Z)-1)]; 
        samplevar=NeweyWest([reshat, mmult( reshat, Z(:,2:cols(Z)) ) ],hac_lag).*scaling;
        if condhomosk==1; samplevar=diag(diag(samplevar)); end; 
        omegahat =  WW*samplevar*WW';
        
      elseif pi>=1;
        lambda=(2/3)*m/n;
        lambdaLZL=1/(2*pi); 
        scaling=[lambda,ones(1,cols(Z)-1).*lambdaLZL;ones(cols(Z)-1,1).*lambdaLZL,ones(cols(Z)-1,cols(Z)-1)]; 
        samplevar=NeweyWest([reshat, mmult( reshat, Z(:,2:cols(Z)) ) ],hac_lag).*scaling;
        if condhomosk==1; samplevar=diag(diag(samplevar)); end; 
        omegahat =  WW*samplevar*WW';
      end; 
      
    % -------------------------------- RECURSIVE SCHEME --------------------
  elseif scheme == 3;
    lambda=1;
    lambdaLZL=(1/pi)*log(1+pi); 
    scaling=[lambda,ones(1,cols(Z)-1).*lambdaLZL;ones(cols(Z)-1,1).*lambdaLZL,ones(cols(Z)-1,cols(Z)-1)]; 
    samplevar=NeweyWest([reshat, mmult( reshat, Z(:,2:cols(Z)) ) ],hac_lag).*scaling;
    if condhomosk==1; samplevar=diag(diag(samplevar)); end; 
    omegahat =  WW*samplevar*WW';
  end;

%==================================================================
% Wald test
teststat = n*(R*deltahat')'*inv(R*omegahat*R')*(R*deltahat');  
pvalue = 1 - cdf('chi2',teststat,rows(R));

%==================================================================
%Prediction interval
i1=Znotdemeaned*deltahat'; 
i2=Znotdemeaned*deltahat'-norminv(1-0.05,0,1)*sqrt(diag(Znotdemeaned*(omegahat/n)*Znotdemeaned'));
interval=[i1,i2];