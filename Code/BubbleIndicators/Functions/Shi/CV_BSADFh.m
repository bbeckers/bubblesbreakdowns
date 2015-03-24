%**************************************************************************
%   "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)
    
%   In this program, we calculate critical value sequences for the backward SADF
%   statistic sequence.
% *************************************************************************

clear all
close all
clc
 
qe=0.95;

m=2499;
T=670; % change your sample size here
r0=0.2;
swindow0=floor(r0*T);    % change your minimum window size here

T00=swindow0:1:T;
dim=size(T00,2);
quantile_badfs=ones(dim,1);

tic
for k=1:1:dim;
    k
    T0=T00(k);
    SI=1;
    randn('seed',SI);   
    e=randn(T0,m); 
    a=T0^(-1);
    y=cumsum(e+a);


    %% THE SUP ADF TEST %%%%%%

    badfs=zeros(T0-swindow0+1,m); 
    sadfs=ones(m,1);

    for j=1:1:m; 
        for i=swindow0:1:T0; 
            badfs(i-swindow0+1,j)= ADF_FL(y(1:i,j),0,1);
        end;  

    sadfs(j,1)=max(badfs(:,j));
    end;

    quantile_badfs(k)=quantile(sadfs,qe);
end;
toc

%% Smooth quantiles over k
% MA (forward and backward) order
q = 2;
quantile_badfs_MA = zeros(dim,1);
for j=1:dim
    quantile_badfs_MA(j) = mean(quantile_badfs(max(j-q,1):min(j+q,dim)));
end
plot(quantile_badfs)
hold on
plot(quantile_badfs_MA,'r')

save(strcat('cv_bsadf_',num2str(T),'_02.mat'))