function pval = CW_test_nan(forec,true,qn);
addpath C:\brossi\research\library\hac;

%This file calculates the p-value of Clark and West's (JoE) test for martingale difference
%hypothesis with the bias correction.

% INPUT: forec is the ROLLING estimate of the forecast and true is realized value; qn is bandwidth in NW
% OUTPUT: pval, the p-value of the test
% This is valid with NaN values
% Test is to distinguish: y(t)=e(t) vs y(t)=x*beta+e(t) where e(t) is mds

a=mean(isfinite(forec),2); c=isfinite(true); b=isfinite(ones(rows(true),1)); 
abc=a+b+c; abc=find(abc==3);
forec=forec(abc,:); true=true(abc); 
kk1=cols(forec); 
teststatv=[]; pval=[]; 
for i1=1:kk1;
        n = length(forec(:,i1));
        y=true.^2-(  (true-forec(:,i1)).^2  - (forec(:,i1)).^2   );
        teststat = sqrt(n)*mean(y)/sqrt(nw(y,qn)); 
        teststatv = [teststatv; teststat]; 
        pval = [pval; 1-cdf('norm',teststat,0,1)];  
end;