function result=Linex(e,a); 
%given a vector of T*1 errors, this function calculates their Linex Loss
%value, where a is the parameter in the Linex Loss 

L=exp(e.*a)-e.*a-ones(rows(e),1);
result=sum(L)/rows(L); 