function [delta,check,AIC,ADF] = ols_adf(y,j,alpha)
%% Header
% This function runs an OLS estimation of an autoregressive model with j
% lagged differences of the level variable y and performs an ADF-test.
%
% Inputs:
% y: Tx1 time series of level variable y
% K: Number of lags to be included of difference of y
% alpha: Significance level to determine the optimal number of lags K by
% subsequent elimination of insignificant lags
%
% Output:
% b: (K+2)x1 vector of estimated parameters
% check: (Kx1) vector with indicators taking the value 1 if the K'th lag is
% insignificant and 0 otherwise.
% AIC: Akaike-Information Criterion for K lags.
% ADF: ADF test statistic

%% Function
if nargin<3
    alpha = 0.05;
end

T = size(y,1);

X = zeros(T-j-1,j);
for i=1:j
    X(:,i) = y(j+1:end-1)-y(j+1-i:end-1-i);
end

if j>0
    X = [ones(T-j-1,1), y(j+1:end-1), X];
else
    X = [ones(T-j-1,1), y(j+1:end-1)];
end
K = j+2;
y = y(j+2:end);

% OLS estimator
b = (X'*X)^(-1)*X'*y;
e = y - X*b;
sigma = e'*e/(T-K-1);
Sigma_b = sigma*(X'*X)^(-1);
t_b = abs(b)./sqrt(diag(Sigma_b));
p_b = 2*(1-tcdf(t_b,T-K-2));

check = p_b(end)<alpha;
AIC = log(sigma)+2*K/T;

ADF = (b(2)-1)./sqrt(Sigma_b(2,2));
delta = b(2);