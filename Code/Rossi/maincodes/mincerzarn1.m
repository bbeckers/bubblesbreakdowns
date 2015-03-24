function result = mincerzarn1(forecerror, yf);

b=olsbeta(forecerror,yf);%[yf,ones(rows(yf),1)]);
Vb=olsvar(forecerror,yf);%[yf,ones(rows(yf),1)]);
wald=b*(inv(Vb))*b';
result=(1 - cdf('chi2',wald,1  )); 