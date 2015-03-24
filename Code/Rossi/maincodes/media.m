function coeff = media(y)

%Calculates mean for a vector or matrix (avg over colmns)
[n k]=size(y);
if n==1 y=y'; [n k]=size(y); end;
somma=sum(y);
coeff=somma./n;



