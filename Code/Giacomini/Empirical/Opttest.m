function result = Opttest(DeltaLt,sigma2,ncountry,model,country_name,tds);
%INPUTS: DeltaLt = sequence of Loss Differences; 
%        sigma2  = an estimate of their variance;
%OUTPUT: the p-value of the One-time test and the p-value of the break test on LM2. 

load pvqlropt.txt; load pvqlrsb.txt;

T=length(DeltaLt); 
t2v=round(T*0.15):1:round(T*0.85); 
LLR7v=[]; LM2v=[]; maxLM2=0; 
   
   for tbreak=t2v;
       LM1=(sum(DeltaLt(1:tbreak)))^2/(T*sigma2);
       LM2=(1/sigma2)*(1/T)*(1/(tbreak/T))*(1/(1-tbreak/T)) ....
           *(sum(DeltaLt(1:tbreak))-(tbreak/T)*sum(DeltaLt(1:T)))^2;
       Chowopt=LM1+LM2;
       LM2v=[LM2v,LM2];
       LLR7v=[LLR7v;Chowopt];
       if LM2>maxLM2; maxLM2=LM2; timebreak=tbreak; end;
   end;
   SupLRopt=max(LLR7v); 
   SupLR=max(LM2v);
   
% P-value
pvSupLRopt=pvcalc(SupLRopt,pvqlropt,1);
pvLM2=pvcalc(SupLR,pvqlrsb,1);
result=[pvSupLRopt;pvLM2]; 
kkk=rows(LM2v'); tplot=[[T*0.15]:1:[T*0.85]]; 
if rows(tplot')==kkk-1; tplot=[tplot,[T*0.85]+1]; 
end;
DeltaLtTrim=DeltaLt(T*0.15:T*0.85,:); time=tds(T*0.15:T*0.85+1,:);
DeltaLt1=mean(DeltaLt(1:timebreak-1,:))*ones(timebreak-[T*0.15],1); 
DeltaLt2=mean(DeltaLt(timebreak:end,:))*ones(cols(tplot)-rows(DeltaLt1),1);% [T*0.85]-timebreak+1,1); 
DeltaLtplot=[DeltaLt1;DeltaLt2];  
if pvLM2<0.10; 
    plot(time,DeltaLtplot);  
    titlestr=['title(''One-break DMW in ',country_name(ncountry,:),' and model ',num2str(model),''');']; eval(titlestr);
    tbreaksaved=time(round(timebreak-[T*0.15]),:); save tbreaksaved tbreaksaved -ascii;
end;