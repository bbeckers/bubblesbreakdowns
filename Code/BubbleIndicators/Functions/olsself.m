function [b,res,BIC,Z] = olsself(y,p)
%% Header
% This function runs an OLS estimation of an autoregressive model with K
%
% Inputs:
% y: Tx1 time series of level variable y
% p: Number of lags to be included of difference of y
%
% Output:
% b: (p+1)x1 vector of estimated parameters
% res: (T-p)x1 vector of residuals
% BIC: Bayesian Information Criterion for p lags.
% Z: (T-p)x(p+1) matrix of regressors

%% Function
T = size(y,1);

Z = zeros(T-p,p);
for i=1:p
    Z(:,i) = y(p+1-i:end-i);
end
Z = [ones(size(Z,1),1),Z];
y = y(p+1:end);

% OLS estimator
b = (Z'*Z)^(-1)*Z'*y;
res = y-Z*b;
sigma_u = res'*res/(T-p-2);
% Information criteria for determining the optimal lag length
BIC = T*log(sigma_u)+(p+1)*log(T);