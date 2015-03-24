function [teststat,pval] = FF_test_c(f_oos_roll,f_in_roll,lambdall,p,m,lambdahh,sigma2,factor,bw);

% Performs the forecast failure test CORRECTED for OVERFITTING
% INPUT: lossdiff, the nx1 sequence of differences between the out-of-sample loss and the in-sample average loss
% lambdall= correction for the losses
% p = # parameters estimated
% m = # in-sample observations
% lambdahh= correction for the h_t
% sigma2=variance of h_t (i.e. the error term in the regression)
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
sigma = sqrt(lambdall*SLL);

teststat = sqrt(n)*mean(lossdiff)/sigma;
teststat = teststat-2*sqrt(n)*(p)*factor/sigma; %this works

pval = 1-normcdf(teststat);
