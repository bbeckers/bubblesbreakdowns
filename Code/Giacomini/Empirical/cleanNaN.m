function result=cleanNaN(x); 
result=[];
for j=1:length(x);
if isnan(x(j,:))==0; result=[result;x(j,:)]; 
end;
end;