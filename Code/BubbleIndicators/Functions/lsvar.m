%% VAR(k,p) LS Estimation.
% y(t)= c+A1*y(t-1)+...+Ap*y(t-p)+e

%% INPUT:
% y [t,k]: time x variable
% p [integer]: # of lags

%% OUTPUT:
% B[k,1+(q-1)+(p*k)]: matrix of estimated VAR coeff.
% A[k,(p*k)]: matrix of estimated VAR coeff. excl. constant and dummies. A=[A1, ..., Ap]
% U[k,t-p]: residuals; 
% SIGMA[k,k]: var-cov matrix; 
% c[k,1]: constants

%% Note

function [B,A,pval,HQ,SC,Q]=lsvar(Y,p,c,h)

if nargin<3
    c = 1;
end
if nargin<4
    h = 4*p;
end

%set up regressors and regressand
[T,K]=size(Y);

y = Y(p+1:T,:)';
if c==1
    X = ones(T-p,1);
else X = [];
end
for i=1:p
    X=[X, Y(p-i+1:T-i,:)];
end;
X = X';

% Run LS regression
B = y*X'*(X*X')^(-1);
vecB = reshape(B,[],1);
if c==1
    A = [B(:,2:end); [eye(K*(p-1)),zeros(K*(p-1),K)]];
else A = [B; [eye(K*(p-1)),zeros(K*(p-1),K)]];
end
U = y-B*X;
SIGMA = U*U'/(T-p-p*K-1);
SIGMA_B = kron((X*X')^(-1),SIGMA);
tstat = reshape(vecB./sqrt(diag(SIGMA_B)),K,[]);
pval = 2*(ones(size(tstat))-tcdf(abs(tstat),T-p-p*K-1));


% Check lag order
SIGMALS = SIGMA*(T-p-p*K-1)/T;
% FPE = ((T+K*p+1)/(T-K*p-1))^K*det(SIGMALS);
% AIC = log(det(SIGMALS))+(2*p*K^2)/T;
HQ = log(det(SIGMALS))+(2*log(log(T))/T)*p*K^2;
SC = log(det(SIGMALS))+(2*log(T)/T)*p*K^2;

% Portmanteau test for residual autocorrelation
U = U';
C = cell(h+1,1);
for j=1:h+1
    C{j,1} = U(j:T-p,:)'*U(1:T-p-j+1,:)/(T-p-j+1);
end
Q = 0;
for j=1:h
    Q = Q+trace(C{j+1,1}'*C{1,1}^(-1)*C{j+1,1}*C{1,1})/(T-p-j);
end
Q = (T-p)^2*Q;
Q = [Q,1-chi2cdf(Q,K^2*(h-p))];