function result = cdfchic(x,k);
%x=asciss of the chisquare
%k=n. degrees of freedom
%This fn calculates the p-values for a chi-square of k degrees of freedom
%It requires table "cdfchictable" in same directory 
%-- that table is taken from Goldberger, Introductory econometrics, 1998, table A.2, p.236
%a copy of the table is in the file "chi_table_list"
%x=7.3;k=5;

load cdfchictable; 
pvalue=[.05:.05:.95,.975,.99,.995]';
if k>13; disp('k too big, increase table'); end;
xvalue=cdfchictable(k,:)';
if x>xvalue(length(xvalue)-1,1); result=0; 
elseif x<=xvalue(1,1); result=1; 
else result=1-cross(pvalue,x*ones(length(pvalue),1),xvalue);
end;

