function res = gmmres(y,z,x,heter);
%This fn. calculates coefficient vector for GMM with optimal weighting matrix (See Hayashi, chp.3)
%y is (T*1) dependent variable
%z is (T*l) independent variable
%x is (T*k) matrix of instruments
%heter is a dummy: =0 if want cond. homosked. and =1 if cond. heterosked.

%e.g. z=randn(100,2);x=z*[2,1;3,4]+randn(100,2);y=z*[1;1]+randn(100,1)*2;heter=0;

%check order condition
T=length(y); Sxx=x'*x/T; Sxz=x'*z/T; Sxy=x'*y/T; 
[a,k]=size(x); [a,l]=size(z);
if k<l; disp('order condition violated'); 
else;
%1st stage (2SLS)
W1=inv(Sxx);
b1=(inv(Sxz'*W1*Sxz))*(Sxz'*W1*Sxy);
e=y-z*b1;
S=zeros(k,k); for i=1:T; S=S+(x(i,:)'*x(i,:)).*(e(i)^2); end; S=S./T;
W2=inv(S);
%2nd stage
b2=(inv(Sxz'*W2*Sxz))*(Sxz'*W2*Sxy);
avarb=inv(Sxz'*W2*Sxz);
if heter==1; b=b2; elseif heter==0; b=b1; end;
end;
res=y-z*b;