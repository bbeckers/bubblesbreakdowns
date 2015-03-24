function res=calds2n(calds,a,b);
uno=cumsum(ones(length(calds),1));
augcalds=[calds,uno]; v=[ones(length(calds),1)*[a,b],uno];
res=intersect(augcalds,v,'rows');
res=res(:,3);