function result = testsoos3(e1,p1,startyear,startmonth,tvec,tds,m,ncountry,model,country_name);
%This file implements the Giacomini-Rossi's One-time test
%m is the size of the window to smooth forecasts

        %clean NaN's
        k=cols(p1); 
        X=[e1,p1,tvec,tds]; X=cleanNaN(X); 
        e1=X(:,1); p1=X(:,2:k+1); time=X(:,end-2:end-1); tds=X(:,end);  %I get NaNs only if all data vector is empty
      
%Oos tests for predictive ability
T=rows(e1); R=calds2n(tvec,startyear,startmonth); 
Pred=T-R;  
rw1=zeros(T,1); true1=zeros(T,1); yoroll1=zeros(T,1);
brollvec=[];  
for j=1:Pred;
   p1here=p1(j:R+j-1,:); 
   broll=olsbeta(e1(j:R+j-1,1),[ones(R,1),p1here]); 
   yoroll1(R+j,1)=[1,p1(R+j,:)]*broll';
   rw1(R+j,1)=0;
   true1(R+j,1)=e1(R+j);
end;
rw=rw1(R+1:T,1);true=true1(R+1:T,1); 
yoroll=yoroll1(R+1:T,1); u1roll=(true-yoroll).^2; u2roll=(true-rw).^2; froll=u2roll-u1roll; 
DeltaLt=u2roll-u1roll; sigma2=cov(DeltaLt); 

tvechere=tvec(R+1:T,:); 
tdshere=tds(R+1:T,:);    mu=m/Pred; 
alpha=0.10; %one-sided

result=Opttest(DeltaLt,sigma2,ncountry,model,country_name,tdshere); 
