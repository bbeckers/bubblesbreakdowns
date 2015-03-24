clear; 
addpath C:\brossi\research\raffaella\march05\MonteCarlo\maincodes; 
addpath C:\brossi\research\raffaella\march05\MonteCarlo; 
delete checkunb.out; diary checkunb.out; 

tableW=[]; tableCONST=[]; tableBETA=[];

for m=[50,100,150,200]; for n=[50,100,150,200]; MC=5000; 
rejWseq=0; rejWsplit=0; rejWroll=0; 
rejtsplit1=0; rejtseq1=0; rejtroll1=0; 
rejtsplit2=0; rejtseq2=0; rejtroll2=0; 
rejWseqC=0; rejWsplitC=0; rejWrollC=0; 
rejtsplit1C=0; rejtseq1C=0; rejtroll1C=0; 
rejtsplit2C=0; rejtseq2C=0; rejtroll2C=0; 

    for rep=1:MC;
        x=[ones(m+n,1),randn(m+n,1)]; y=randn(m+n,1); 
        result=linear_FF1e_new(y,x,m,1,'forc');
        true=result(:,4); FFseq=true-result(:,1); FFsplit=true-result(:,2); FFroll=true-result(:,3);
        z=randn(n,1); regr=[ones(n,1),z]; 
        UNBthetaSEQ=olsbeta(FFseq,regr)'; 
        UNBthetaSPLIT=olsbeta(FFsplit,regr)';
        UNBthetaROLL=olsbeta(FFroll,regr)';
        VthetaSEQ=cov(olsres(FFseq,regr)')*inv(regr'*regr/n); 
        VthetaSPLIT=cov(olsres(FFsplit,regr)')*inv(regr'*regr/n);
        VthetaROLL=cov(olsres(FFroll,regr)')*inv(regr'*regr/n);
        Wseq=n*UNBthetaSEQ'*(inv(VthetaSEQ))*UNBthetaSEQ; 
        Wsplit=n*UNBthetaSPLIT'*(inv(VthetaSPLIT))*UNBthetaSPLIT; 
        Wroll=n*UNBthetaROLL'*(inv(VthetaROLL))*UNBthetaROLL; 
        tseq1=sqrt(n)*UNBthetaSEQ(1,1)'*(inv(sqrt(VthetaSEQ(1,1))));  
        tsplit1=sqrt(n)*UNBthetaSPLIT(1,1)'*(inv(sqrt(VthetaSPLIT(1,1)))); 
        troll1=sqrt(n)*UNBthetaROLL(1,1)'*(inv(sqrt(VthetaROLL(1,1)))); 
        if cols(regr)>1; 
        tseq2=sqrt(n)*UNBthetaSEQ(2,1)'*(inv(sqrt(VthetaSEQ(2,2)))); 
        tsplit2=sqrt(n)*UNBthetaSPLIT(2,1)'*(inv(sqrt(VthetaSPLIT(2,2)))); 
        troll2=sqrt(n)*UNBthetaROLL(2,1)'*(inv(sqrt(VthetaROLL(2,2)))); 
        end;
    
        lambda=1+n/m; pai=n/m;
        if Wseq>5.99; rejWseq=rejWseq+1; end; 
        if Wsplit>5.99; rejWsplit=rejWsplit+1; end; 
        if Wroll>5.99; rejWroll=rejWroll+1; end; 
        if abs(tseq1)>1.96; rejtseq1=rejtseq1+1; end; 
        if abs(tsplit1)>1.96; rejtsplit1=rejtsplit1+1; end; 
        if abs(troll1)>1.96; rejtroll1=rejtroll1+1; end; 
        if cols(regr)>1;
            if abs(tseq2)>1.96; rejtseq2=rejtseq2+1; end; 
            if abs(tsplit2)>1.96; rejtsplit2=rejtsplit2+1; end; 
            if abs(troll2)>1.96; rejtroll2=rejtroll2+1; end; 
        end; 
        %CORRECTED TEST
        SS=cov(olsres(FFseq,regr)'); sigma=SS(1,1); D=zeros(cols(x)-1,1); B=x(:,2)'*x(:,2)/rows(x); h=mmult(olsres(y(1:m),x(1:m,:)),x(1:m,2)); l=FFseq; HAC_lag=0;
        VthetaSEQ = varPE(z,sigma,D,B,h,3,l,HAC_lag,pai);
        SS=cov(olsres(FFsplit,regr)'); sigma=SS(1,1); D=zeros(cols(x)-1,1); B=x(:,2)'*x(:,2)/rows(x); h=mmult(olsres(y(1:m),x(1:m,:)),x(1:m,2)); l=FFsplit; HAC_lag=0;
        VthetaSPLIT = varPE(z,sigma,D,B,h,1,l,HAC_lag,pai);
        SS=cov(olsres(FFroll,regr)'); sigma=SS(1,1); D=zeros(cols(x)-1,1); B=x(:,2)'*x(:,2)/rows(x); h=mmult(olsres(y(1:m),x(1:m,:)),x(1:m,2)); l=FFroll; HAC_lag=0;
        VthetaROLL = varPE(z,sigma,D,B,h,2,l,HAC_lag,pai);
        WseqC=n*UNBthetaSEQ'*(inv(VthetaSEQ))*UNBthetaSEQ; 
        WsplitC=n*UNBthetaSPLIT'*(inv(VthetaSPLIT))*UNBthetaSPLIT; 
        WrollC=n*UNBthetaROLL'*(inv(VthetaROLL))*UNBthetaROLL; 
        tseq1C=sqrt(n)*UNBthetaSEQ(1,1)'*(inv(sqrt(VthetaSEQ(1,1))));  
        tsplit1C=sqrt(n)*UNBthetaSPLIT(1,1)'*(inv(sqrt(VthetaSPLIT(1,1)))); 
        troll1C=sqrt(n)*UNBthetaROLL(1,1)'*(inv(sqrt(VthetaROLL(1,1)))); 
        if cols(regr)>1; 
        tseq2C=sqrt(n)*UNBthetaSEQ(2,1)'*(inv(sqrt(VthetaSEQ(2,2)))); 
        tsplit2C=sqrt(n)*UNBthetaSPLIT(2,1)'*(inv(sqrt(VthetaSPLIT(2,2)))); 
        troll2C=sqrt(n)*UNBthetaROLL(2,1)'*(inv(sqrt(VthetaROLL(2,2)))); 
        end;
    
        lambda=1+n/m; 
        if WseqC>5.99; rejWseqC=rejWseqC+1; end; 
        if WsplitC>5.99; rejWsplitC=rejWsplitC+1; end; 
        if WrollC>5.99; rejWrollC=rejWrollC+1; end; 
        if abs(tseq1C)>1.96; rejtseq1C=rejtseq1C+1; end; 
        if abs(tsplit1C)>1.96; rejtsplit1C=rejtsplit1C+1; end; 
        if abs(troll1C)>1.96; rejtroll1C=rejtroll1C+1; end; 
        if cols(regr)>1;
            if abs(tseq2C)>1.96; rejtseq2C=rejtseq2C+1; end; 
            if abs(tsplit2C)>1.96; rejtsplit2C=rejtsplit2C+1; end; 
            if abs(troll2C)>1.96; rejtroll2C=rejtroll2C+1; end; 
        end;

        
    end;
    tableW=[tableW; m,n,rejWseq/MC,rejWsplit/MC,rejWroll/MC,rejWseqC/MC,rejWsplitC/MC,rejWrollC/MC];
    tableCONST=[tableCONST; m,n,rejtseq1/MC, rejtsplit1/MC, rejtroll1/MC, rejtseq1C/MC, rejtsplit1C/MC, rejtroll1C/MC],
    tableBETA=[tableBETA; m,n,rejtseq2/MC, rejtsplit2/MC, rejtroll2/MC, rejtseq2C/MC, rejtsplit2C/MC, rejtroll2C/MC];
    
    %disp([m,n,rejWseq/MC,rejWsplit/MC,rejWroll/MC, rejtseq1/MC, rejtsplit1/MC, rejtroll1/MC, rejtseq2/MC, rejtsplit2/MC, rejtroll2/MC]);
    %disp([m,n,rejWseqC/MC,rejWsplitC/MC,rejWrollC/MC, rejtseq1C/MC, rejtsplit1C/MC, rejtroll1C/MC, rejtseq2C/MC, rejtsplit2C/MC, rejtroll2C/MC]);
end; end;  

disp(tableCONST);
disp(tableBETA);
disp(tableW);
diary off;