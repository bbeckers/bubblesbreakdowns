function [teststat,pval] = FF_test(f_oos_roll,f_in_roll,lambdahh,bw);

% Performs the forecast failure test
% INPUT: lossdiff, the nx1 sequence of differences between the out-of-sample loss and the in-sample average loss
% OUTPUT: pval, the p-value of the forecast failure test (rejection indicates forecast failure)   
%
% ASSUMPTIONS: lossdiff is iid

lossdiff=f_oos_roll-f_in_roll; 
n = length(lossdiff);

LOtilde = lossdiff-mean(lossdiff);

% HAC variance estimator of demeaned surprise losses
if bw==0
    SLL = cov(lossdiff);
else
    SLL = 0;
    for j=1:bw
        Gamma = LOtilde(1+j:end)'*LOtilde(1:end-j)/n;
        SLL = SLL+2*(1-j/(bw+1))*Gamma;
    end
    SLL = SLL+LOtilde'*LOtilde/n;
end
% Variance estimator
sigma = sqrt(lambdahh*SLL);

teststat = sqrt(n)*mean(lossdiff)/sigma;
pval = 1-normcdf(teststat);
