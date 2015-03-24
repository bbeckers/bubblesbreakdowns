function pvalue=pvcalc(osser,tavola,kk);
%Calculates p-values of TVP and Optimal tests -- requires tables in the same directory
%Inputs: osser=osserved value of the statistic; tavola=table with critical values; kk=number of regressors
%Output: p-value (scalar)

numero=tavola(:,[1,kk+1+1]);
uno=0; 
if osser<=tavola(1,kk+2) pv=1; uno=uno+1; end;
if osser>=tavola(34,kk+2) pv=0; uno=uno+1; end;
if uno==0;
   rigal=find(numero(:,2)<=osser);
   riga=[max(rigal);max(rigal)+1]; sel=tavola(riga,[1,kk+1+1]);
   pv=sel(2,1)+(sel(2,2)-osser)*(sel(1,1)-sel(2,1))/(sel(2,2)-sel(1,2));
   end;
pvalue=pv;
