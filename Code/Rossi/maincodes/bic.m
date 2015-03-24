function result=bic(y,x,pmax,const);
%Runs ols regression of y onto x and selects p by BIC, where p is the lag
%length of every regressor in x

T=rows(y); SSR=[]; 
for i=pmax:-1:1;
    if const==0; xhere=x(i:T,:); elseif const==1; xhere=[ones(T-i+1,1),x(i:T,:)]; end;
    for k=1:i-1; xhere=[xhere,x(i-k:T-k,:)]; end;
    e=olsres(y(i:T),xhere); SSR=[SSR;i,log((e'*e)/rows(e))+cols(xhere)*(log(rows(e)))/rows(e)]; 
    %i, xhere, SSR
end;
SSRsorted=sortrows(SSR,2); 
result=SSRsorted(1,1); 