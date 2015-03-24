function result = testsoos2(e1,p1,startyear,startmonth,tvec,tds);
%This file calculates Clark-West's (2005) and the Diebold-Mariano-West's 
%traditional tests and p-values

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

%Rolling tests
yoroll=yoroll1(R+1:T,1); u1roll=(true-yoroll).^2; u2roll=(true-rw).^2; froll=u1roll-u2roll; 
varfroll=cov(froll)/(Pred); 
DMroll=media(froll)/sqrt(varfroll);  
%IN PAPELL-MOLODTSOVA THE DMW TEST IS ONE-SIDED:
pvDMroll=1-normcdf(-DMroll,0,1); 

CWroll=Pred*(media(u2roll.^2-u1roll.*u2roll))/(media((u1roll-media(u1roll)).^2)); %disp('ENCNEW roll');disp(oosvcmroll);disp('10%,5%,1% cv');disp([ooscv(kk+1,0.1,2,2),ooscv(kk+1,0.05,2,2),ooscv(kk+1,0.01,2,2)]);

pvCWroll=CW_test_nan1(yoroll,true,0);  

result=[DMroll,CWroll; pvDMroll,pvCWroll];
