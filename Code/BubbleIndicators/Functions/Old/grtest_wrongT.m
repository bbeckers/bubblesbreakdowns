function [SL,SLbar,tc,pc,t,p] = grtest_wrongT(LO,LI,Y,X,scheme,bw)

%% Description %%
% ----------------------------------------------------------------------- %
% This function conducts the forecast breakdown test by Giacomini & Rossi
% (2009, The Review of Economic Studies).
% ----------------------------------------------------------------------- %
% INPUT to the function:
% LO:  [n] Vector of out-of-sample losses
% LI:  [n] Vector of average in-sample losses for each estimation point t
% Y:   [T] Vector of target variable
% X:   [(T)x(k)] Matrix of predictors
% tau: Scalar denoting forecast horizon
% s:   String sececting the forecasting scheme
%      'fixed' (default), 'rolling', 'recursive'
% ----------------------------------------------------------------------- %
% OUTPUT of the function:
% SL:    [n] Vector of surprise losses
% SLbar: Scalar average surprise loss
% tc:    Overfitting-corrected test-statistic
% pc:    Overfitting-corrected p-value
% t:     Uncorrected test statistic
% p:     Uncorrected p-value
% ----------------------------------------------------------------------- %

%% Check inputs for consistency
if size(LO,1)~=size(LI,1)
    error('Dimension of output losses must correspond to estimation points')
end

if nargin<6
    scheme = 'recursive';
    warning('Scheme is unspecified and was set to "recursive" (default)')
elseif ~strcmp(scheme,'fixed') && ~strcmp(scheme,'rolling') && ~strcmp(scheme,'recursive')
    error('Variable s must contain either "fixed", "rolling" or "recursive"')
end

%% Begin function
n = size(LO,1);
[T,k] = size(X);
m = T-n;

% Adjustment paramater for surprise loss variance estimator
if strcmp(scheme,'fixed')
    lambda = 1+n/m;
elseif strcmp(scheme,'rolling') && n<m
    lambda = 1-(1/3)*(n/m)^2;
elseif strcmp(scheme,'rolling') && n>=m
    lambda = (2/3)*(m/n);
else
    lambda = 1;
end

% Adjustment paramater in overfitting correction
if strcmp(scheme,'recursive')
    gamma = (1/sqrt(n))*log(1+n/m);
else
    gamma = sqrt(n)/m;
end

%% Average surprise losses
% Vector of surprise losses
SL = LO-LI;
% Out-of-sample mean of surprise losses
SLbar = mean(SL);

%% Asymptotic variance
% Demeaned surprise losses
LOtilde = LO-mean(LO);
% HAC variance estimator of demeaned surprise losses
if bw==0
    SLL = cov(LO);
else
    SLL = 0;
    for j=1:bw
        Gamma = LOtilde(1+j:end)'*LOtilde(1:end-j)/n;
        SLL = SLL+2*(1-j/(bw+1))*Gamma;
    end
    SLL = SLL+LOtilde'*LOtilde/n;
end
% Variance estimator
sigma = sqrt(lambda*SLL);

%% Test statistic
t = sqrt(n)*SLbar/sigma;
p = 1-normcdf(t);

%% Overfitting correction
% Asymptotic variance covariance matrix of coefficients
XX = X'*X;
beta = XX^(-1)*X'*Y;
u = Y-X*beta;
Vbeta = u'*u/(T-k)*XX^(-1);
% Correction parameter
c = 2*gamma*trace(XX/T*Vbeta);
% Adjusted test statistic
tc = t-c/sigma;
% tc = t-c/sigma;
pc = 1-normcdf(tc);