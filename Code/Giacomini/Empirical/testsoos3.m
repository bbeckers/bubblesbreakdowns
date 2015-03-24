function result = testsoos3(e1,p1,startyear,startmonth,tvec,tds,m);
%This file calculates Giacomini-Rossi's Fluctuation test with both the Clark-West's (2005) and the Diebold-Mariano-West's 
%statistics, and plots the Giacomini-Rossi's Fluctuations bands
 
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
%Rolling tests
DMrollv=[]; CWrollv=[]; 
for s=1:Pred-m+1;
    frollhere=froll(s:s+m-1,:); 
    varfrollhere=cov(frollhere)/m; 
    DMrollhere=media(frollhere)/sqrt(varfrollhere);  
    DMrollv=[DMrollv;DMrollhere]; 
    CWroll=CW_test_nanJAE(yoroll(s:s+m-1,:),true(s:s+m-1,:),0);  
    CWrollv=[CWrollv;CWroll];
end; 

tvechere=tvec(R+1+round(m/2):end,:); tvechere=tvechere(1:Pred-m+1,:); 
tdshere=tds(R+1+round(m/2):end,:);   tdshere=tdshere(1:Pred-m+1,:); mu=m/Pred; 
alpha=0.10; %one-sided
k=GiacominiRossiCV(mu,alpha);


result=[DMrollv,CWrollv,tvechere,tdshere,k*ones(rows(DMrollv),1),-k*ones(rows(DMrollv),1)];
