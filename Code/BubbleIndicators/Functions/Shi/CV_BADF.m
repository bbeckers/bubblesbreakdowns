%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)
    
%   In this program, we calculate critical value sequences for the backward ADF
%   statistic sequence.
% *************************************************************************

 

clear all
close all
clc

format short
 
qe=[0.90;0.95;0.99];

tic;

m=5000;
T=1680;            % change your sample size here
swindow0=36;       % change your minimum window size here

T00=swindow0:1:T;
dim=size(T00,2);
quantile_badfs=ones(3,dim);

for k=1:1:dim;
   T0=T00(k);

   SI=1;
   randn('seed',SI);   
   e=randn(T0,m); 
   a=T0^(-1);
   y=cumsum(e+a);

  adfs=ones(m,1);
  for j=1:1:m; 
      adfs(j)= ADF_FL(y(:,j),0,1);  
  end;

quantile_badfs(:,k)=quantile(adfs,qe);
end;

dlmwrite('CV_BADF1680_36.txt', [qe quantile_badfs],'-append', 'delimiter', '\t','precision', '%14.7f')
type CV_BADF1680_36.txt
