function result = olsstddevnw(y,x,qn)

%Calculates residuals for an unrestricted VAR(p)
%y must have different variables as columns
%variance is corrected for small samples (n-k) (18-3-03)

[n k]=size(y); [nx,kx]=size(x);
phi=zeros(k,1); res=zeros(n,1);
coeff=((inv(x'*x))*(x'*y))';
res=y-x*coeff';
sigma2=nw(res,qn);
%sigma2=(res'*res)/(length(res)-kx);
length(res)-kx;
var=(inv(x'*x))*sigma2;
result=sqrt(diag(var));
