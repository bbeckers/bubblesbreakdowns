function [teststat,pval] = FF_test_McC3(f_oos_roll,f_in_roll,lambdahh);

% Performs the forecast failure test
% INPUT: lossdiff, the nx1 sequence of differences between the out-of-sample loss and the in-sample average loss
% OUTPUT: pval, the p-value of the forecast failure test (rejection indicates forecast failure)   
%
% ASSUMPTIONS: lossdiff is iid

lossdiff=f_oos_roll-f_in_roll; 
n = length(lossdiff);
teststat = sqrt(n)*mean(lossdiff)/sqrt(cov([f_oos_roll])*lambdahh);
pval = 2*(1 - cdf('norm',abs(teststat),0,1));
