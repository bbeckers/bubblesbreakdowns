function res = lagnmatrix(x,p);
%Lagging a vector x p times --- missing values are cut
%  say x=[1 2 3 4 5];
% p=2 you get 
% 0 0 
% 1 0
% 2 1
% 3 2
% 4 3
 
T=rows(x);
res=[zeros(p,1);x(1:T-p,1)]; 
for j=1:p-1; res=[[zeros(p-j,1);x(1:T-p+j,1)],res]; end; %res=[[zeros(p-j,1);x(1:T-j,1)],res]; end;
%res=res(p+1:T,:);