function result=rollavg(y,q);
%calculates rolling future averages of y with a q window, leaving the last q-1
%values unchanged;
result=[];
for j=1:rows(y)-q+1;
    result=[result;sum(y(j:q+j-1,:))/q];
end;
result= [ result;y(rows(y)-q+2:rows(y),:)];     