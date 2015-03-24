function result=calendarm(fyds,fmds,lyds,lmds);

nds =   12*(lyds-fyds)  + (lmds-fmds)+1;  % Data Set Sample Size @
calds=zeros(nds,2); calds(1,1)=fyds; calds(1,2)=fmds;
yr=fyds; mt=fmds;
i=2; for i=2:nds;
 mt=mt+1;
 if mt > 12; mt=1; yr=yr+1; end; 
 calds(i,1)=yr; calds(i,2)=mt;
 end; 
 result=calds;
