function coeff = olsbeta(y,x)

%Calculates ols estimates

[n k]=size(y);
phi=zeros(k,1); res=zeros(n,1);
coeff=((inv(x'*x))*(x'*y))';

