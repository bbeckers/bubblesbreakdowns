%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)
    
%   In this program, we calculate critical values for the sup 
%   ADF statistic.
% *************************************************************************

 

clear all
close all
clc

format short
 
qe=[0.90;0.95;0.99];

tic;

m=5000;
T=1680;          % change your sample size here
swindow0=36;     % change your minimum window size here
r0=36/T;  

%% %%%% DATA GENERATING PROCESS %%%%%%
SI=1;
randn('seed',SI);   
e=randn(T,m); 
a=T^(-1);
y=cumsum(e+a);

%% THE SUP ADF TEST %%%%%%

badfs=zeros(T-swindow0+1,m); 
sadf=ones(m,1);
 
for j=1:1:m; 
    for i=swindow0:1:T; 
      badfs(i-swindow0+1,j)= ADF_FL(y(1:i,j),0,1);  
    end;   
  end;
  
sadf(:,1)=max(badfs,[],2);


quantile_sadf=quantile(sadf,qe);

dlmwrite('CV_BADF1680_36.txt', quantile_sadf,'-append', 'delimiter', '\t','precision', '%14.7f')
type CV_BADF1680_36.txt
