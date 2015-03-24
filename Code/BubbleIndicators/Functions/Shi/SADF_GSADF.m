%*********************************************************************
%    "Testing for Multiple Bubbles" by Phillips, Shi and Yu (2011)
    
%    In this program, we calculate the sup ADF statistic 
%    (the backward ADF statistic sequence) and the generalized 
%    sup ADF statistic (the backward SADF statistic sequence) 
% *******************************************************************


clear all
close all
clc

format short
 
SPDV=xlsread('C:\Research\SP_DV.xlsx');  
y=SPDV;

T=length(y);
swindow0=36;
r0=36/T;  
dim=T-swindow0+1;


%% THE SUP ADF TEST %%%%%%

badfs=zeros(T-swindow0+1,1); 
for i=swindow0:1:T; 
  badfs(i-swindow0+1,1)= ADF_FL(y(1:i,1),0,1);  
end;   
sadf=max(badfs);

display('The sup ADF statistic is'); sadf

figure (1)
plot(badfs);
title('The backward ADF sequence','FontSize',10);

pause

 %% THE GENERALIZED SUP ADF TEST %%%%%%
r2=swindow0:1:T;
r2=r2';
rw=r2-swindow0+1;

bsadfs=zeros(1,dim); 

for v=1:1:size(r2,1);
    swindow=swindow0:1:r2(v);
    swindow=swindow';
    r1=r2(v)-swindow+1;
    rwadft=zeros(size(swindow,1),1);
    for i=1:1:size(swindow,1); 
       rwadft(i)= ADF_FL(y(r1(i):r2(v),1),0,1);   % two tail 5% significant level
    end;  
    bsadfs(1,v)=max(rwadft);
 end;
  
gsadf=max(bsadfs(1,:)');



display('The generalized sup ADF statistic is'); gsadf

figure (2)
plot(bsadfs);
title('The backward SADF sequence','FontSize',10);
 

