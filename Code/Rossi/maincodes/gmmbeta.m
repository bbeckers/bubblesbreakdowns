function b = gmmbeta(y,z,x,heter);
%Calculates coefficient vector for GMM with optimal weighting matrix (See Hayashi, chp.3)
%y is (T*1) dependent variable
%z is (T*l) independent variable
%x is (T*k) matrix of instruments
%heter is a dummy: =0 if want cond. homosked. and =1 if cond. heterosked.
%we checked that if x=z then we get back to the OLS if cond homosk
if isempty(x)==1; x=[]; end; if isempty(z)==1; z=[]; end; 
%check order condition
T=length(y); Sxx=x'*x/T; Sxz=x'*z/T; Sxy=x'*y/T; [a,k]=size(x); [a,l]=size(z);
if k<l; disp('order condition violated'); 
else;
%1st stage (2SLS)
W1=inv(Sxx);
b1=(inv(Sxz'*W1*Sxz))*(Sxz'*W1*Sxy);
e=y-z*b1;
S=x'*(diag(e.^2))*x./T;
W2=inv(S);
%2nd stage
b2=(inv(Sxz'*W2*Sxz))*(Sxz'*W2*Sxy); b=b2;
end;