function avarb = gmmvar(y,z,x,heter);
%Calculates variance of the coefficient vector for GMM with optimal weighting matrix (See Hayashi, chp.3)
%y is (T*1) dependent variable
%z is (T*l) independent variable
%x is (T*k) matrix of instruments
%heter is a dummy: =0 if want cond. homosked. and =1 if cond. heterosked.
%we checked that if x=z then we get back to the OLS if cond homosk

%check order condition
T=length(y); Sxx=x'*x/T; Sxz=x'*z/T; Sxy=x'*y/T; [a,k]=size(x); [a,l]=size(z);
if k<l; disp('order condition violated'); 
    else;
    %1st stage (2SLS)
    W1=inv(Sxx);
    b1=(inv(Sxz'*W1*Sxz))*(Sxz'*W1*Sxy);
    e=y-z*b1;
    if heter==1; S=x'*(diag(e.^2))*x./T; W2=inv(S);
    elseif heter==0; S=inv(Sxx)*sum(e.^2)/T; W2=inv(S); end;
    %2nd stage
    b2=(inv(Sxz'*W2*Sxz))*(Sxz'*W2*Sxy);
    if heter==1; avarb=inv(Sxz'*W2*Sxz); elseif heter==0; avarb=inv(Sxx)*sum((y-z*b2).^2)/T; end;
end;