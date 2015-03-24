function [Wstat,pvalW,SLfitCI] = FB_Wald_CI(Z,reshat,deltahat,n,m,sigma2,scheme,bw,alpha)

%% Description %%
% ----------------------------------------------------------------------- %
% This function conducts the Wald test on predicting surprise losses by a
% set of explanatory variables and analysing their behavior over time by 
% confidence intervals as proposed by Giacomini & Rossi (2009, The Review
% of Economic Studies).
% ----------------------------------------------------------------------- %
% INPUT to the function:
% Z:        [(tZ)x(k)] Matrix of predictors
% reshat:   [tZ] Vector of regression residuals
% deltahat: [k] vector of paramater estimates
% n:        Evaluation window size
% m:        Training sample length for estimation
% sigma2:   Variance of out-of-sample surprise losses
% scheme:   String sececting the forecasting scheme
%           'fixed' (default), 'rolling', 'recursive'
% bw:       Bandwidth for HAC estimator
% alpha:    Confidence level for CI (alpha=0.05)
% ----------------------------------------------------------------------- %
% OUTPUT of the function:
% Wstat:   Wald test statistic
% pvalW:   p-value for Wald test statistic
% SLfitCI: [tZ] vector providing the CI for fitted surprise losses
% ----------------------------------------------------------------------- %

% Scaling parameter of covariance matrix for Wald Test
pie = n/m;
if strcmp(scheme,'recursive')
    LAMBDA = log(1+pie)/pie;
elseif strcmp(scheme,'rolling') && n<=m
    LAMBDA = 1-pie/2;
elseif strcmp(scheme,'rolling') && n>m
    LAMBDA = (2*pie)^(-1);
else
    LAMBDA = 1;
end

% Demeaned regressors
Ztilde = Z(:,2:end);
[tZ,k] = size(Ztilde);
Ztilde = Ztilde-ones(size(Z,1),1)*mean(Ztilde);

% Inputs to covariance matrix for Wald Test
Szz = zeros(k);
SzLL = zeros(k,1);
SzLzL = zeros(k);

for j=1:bw
%     GAMMAzz = Ztilde(1+j:end,:)'*Ztilde(1:end-j,:)/tZ;
%     Szz = Szz + (1-j/(bw+1))*(GAMMAzz+GAMMAzz');
    GAMMAzLL = (Ztilde(1+j:end,:)'*(reshat(1+j:end).*reshat(1:end-j))+Ztilde(1:end-j,:)'*(reshat(1:end-j).*reshat(1+j:end)))/tZ;
    SzLL = SzLL + (1-j/(bw+1))*GAMMAzLL;
    GAMMAzLzL = (Ztilde(1+j:end,:).*(reshat(1+j:end)*ones(1,k)))'*(reshat(1:end-j)*ones(1,k).*Ztilde(1:end-j,:))/tZ;
    SzLzL = SzLzL + (1-j/(bw+1))*(GAMMAzLzL+GAMMAzLzL');
end
Szz = Szz+Ztilde'*Ztilde/tZ;
SzLL = SzLL+Ztilde'*reshat.^2/tZ;
SzLzL = SzLzL+(Ztilde.*(reshat*ones(1,k)))'*(Ztilde.*(reshat*ones(1,k)))/tZ;

A = [1, -mean(Z(:,2:end))*Szz^(-1); zeros(k,1), Szz^(-1)];
C = [sigma2, LAMBDA*SzLL'; LAMBDA*SzLL, SzLzL];

% Covariance matrix for Wald Test
Omegahat = A*C*A';

% Test statistic and p-value
Wstat = n*deltahat'*Omegahat^(-1)*deltahat;
pvalW = 1-chi2cdf(Wstat,k+1);

% Confidence interval for fitted surprise losses
SLfitCI = Z*deltahat-norminv(1-alpha)*sqrt(diag((Z*Omegahat*Z'/n)));