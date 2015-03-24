function DynPT = pt_test(Y,x,pmax)

% This function performs the Test for Dependence among serially correlated
% multicategory variables by Pesaran and Timmermann (2009)
% Inputs
% Y is a vector of forecast values
% X is a vector of realizations
% pmax is the maximum lag order for the dynamic regression model
% Outputs
% DynPT is the test statistic (column 1) and the corresponding p-value

AIC = zeros(pmax,1);
T = size(Y,1);

% Find most parsimonous, best-fitting model
for p=1:pmax
    X = zeros(T-p,2*p);
    for i=1:p
        X(:,i) = x(p+1-i:end-i);
    end
    for j=1:p
        X(:,p+j) = Y(p+1-j:end-j);
    end
    X = [ones(T-p,1) x((p+1:end)) X];
    y = Y(p+1:end);
    B = (X'*X)^(-1)*X'*y;
    E = y-X*B;
    AIC(p) = (T-p)*log(E'*E/(T-p))+2*p;
end; clear p
[~,pstar] = min(AIC);

% Find canonical correlations for best fitting model
X = zeros(T-pstar,pstar);
for i=1:pstar
    X(:,i) = x(pstar+1-i:end-i);
end; clear i
for j=1:pstar
    X(:,pstar+j) = Y(pstar+1-j:end-j);
end; clear j
% Auxiliary matrices
W = [ones(T-pstar,1) x((pstar+1:end)) X];
Y = Y(pstar+1:end);
M = eye((T-pstar))-W*(W'*W)^(-1)*W';
% Correlations
Syx = Y'*M*X/(T-pstar);
Sxy = X'*M*Y/(T-pstar);
Sxx = X'*M*X/(T-pstar);
Syy = Y'*M*Y/(T-pstar);
% Test statistic
S = (T-pstar-1)*Syy^(-1)*Syx*Sxx^(-1)*Sxy;
DynPT = [S, 1-chi2cdf(S,1)];
