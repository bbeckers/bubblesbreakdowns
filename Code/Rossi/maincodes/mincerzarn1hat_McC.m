function result = mincerzarn1hat_McC(forecerror, yf,lambdahh);

b=olsbeta(forecerror,yf);%[yf,ones(rows(yf),1)]);
Vb=olsvar(forecerror,yf);%[yf,ones(rows(yf),1)]);
wald=b*(inv(Vb*lambdahh))*b';
result=(1 - cdf('chi2',wald,1  )); 