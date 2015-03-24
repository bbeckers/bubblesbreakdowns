%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)
    
%   In this program, we calculate critical values for the generalized sup 
%   ADF statistic.
% *************************************************************************
 

% clear all
% close all
% clc

format short
 
qe=[0.90;0.95;0.99];

tic;

m=100;
T=400;              % change your sample size here
swindow0=40;         % change your minimum window size here
r0=0.1;  
dim=T-swindow0+1;

%% %%%% DATA GENERATING PROCESS %%%%%%
SI=1;
randn('seed',SI);   
e=randn(T,m); 
a=T^(-1);
y=cumsum(e+a);


%% THE GENERALIZED SUP ADF TEST %%%%%%
r2=swindow0:1:T;
r2=r2';
rw=r2-swindow0+1;

gsadf=ones(m,1);
sadfs2=zeros(m,dim); 
  
for j=1:1:m; 

 for r2=swindow0:1:T;
    dim0=r2-swindow0+1;
    rwadft=zeros(dim0,1);
    for r1=1:1:dim0; 
       rwadft(r1)= ADF_FL(y(r1:r2,j),0,1);  % two tail 5% significant level
    end;
    sadfs2(j,r2)=max(rwadft);
 end;
  
gsadf(j,1)=max(sadfs2(j,:),[],2);
j
end;

toc;

quantile_gsadf=quantile(gsadf,qe);


dlmwrite('CV_BSADF1680_36.txt', quantile_gsadf,'-append', 'delimiter', '\t','precision', '%14.7f')

type CV_BSADF1680_36.txt




%
% if simulation==1
%     SI = 1;
%     randn('seed',SI);
%     % Generate random numbers for random walk with drift
%     e = randn(T,sims);
%     a = T^(-1);
%     y = cumsum(e+a);
%     % Start simulation
%     ADFSIM = zeros(sims,dim);
%     for m=1:sims; 
%         for rend=tau0:1:T;
%             dim0 = rend-tau0+1;
%             rwadft = zeros(dim0,1);
%             for r1=1:1:dim0;
%                 [~,~,~,rwadft(r1)] = ols_adf(y(r1:rend,m),0);
%             end;  
%             ADFSIM(m,rend-tau0+1)=max(rwadft);
%         end;
%     end
%     cv_ADF = quantile(ADFSIM,alpha);
% end
