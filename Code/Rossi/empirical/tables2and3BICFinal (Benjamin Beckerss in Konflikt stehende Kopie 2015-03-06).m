
%----------------------------------------------------
%BIC
%----------------------------------------------------
delete empirical15bic.out; diary empirical15bic.out; 

addpath('..\maincodes')
disp('BIC RESULTS');

clear; tableh=[]; ForecBreakd=[]; ForecBreakd1=[]; wanttoplot=0;
hvec=[1,3,12]; Table3new=[]; table=[];
for kk=1:cols(hvec); h=hvec(1,kk),
    yearmonthRv=[1979,1];%1985,1;1993,1];[1984,12; 1978,12; 1976,12; 1975,12; 1993,1]; 

    for yearmonthR=1:rows(yearmonthRv);
    yearR=yearmonthRv(yearmonthR,1); monthR=yearmonthRv(yearmonthR,2); 


    cd data; cd qvmd; load ruc_latest.txt; load cpi_latest.txt; cd .. ; cd .. ;
    cd data; load fedfund5901.txt; load m2.txt; load commodp.txt; load tbill3m.txt; load tb10years.txt; cd .. ;

    y=cpi_latest; y_sum=[zeros(11,1);[100*log(y(13:rows(y),1)./y(1:rows(y)-12,1))]]; 
    y=1200*log(y(2:rows(y),:)./y(1:rows(y)-1,:)); y_init=y; 
    pi=commodp; pi=[0;1200*log(pi(2:rows(pi),:)./pi(1:rows(pi)-1,:))]; 
    m2=[0;log(m2(2:rows(m2),:)./m2(1:rows(m2)-1,:))];
    spread=tbill3m-tb10years; spread=[0;spread(2:rows(spread),:)];
    m=ruc_latest; m=m(2:rows(m),:); 
    %results=ols(ruc_latest,[ones(rows(ruc_latest),1),cumsum(ones(rows(ruc_latest),1)),cumsum(ones(rows(ruc_latest),1)).^2]),
    %gap=results.resid; %gap=gap(2:rows(m),:); 
    gap=-olsres(ruc_latest,[ones(rows(ruc_latest),1),cumsum(ones(rows(ruc_latest),1)),cumsum(ones(rows(ruc_latest),1)).^2]);
    gap=rollavg1(gap,12); gap=gap(2:rows(gap),:); ffr=fedfund5901;  

    msef=0; mseb=0; msegc=0; msewald=0;

    % -- Construct Calendar Sequences for Plotting, etc -- @
    fyds=1948; fmds=2;                      % First Year (fyds) and month (fmds) of Data Set and vintage (fyvt)%
    lyds=2004; lmds=6;                      % Last Year of Data Set and last Month of Data Set %

    nds =   12*(lyds-fyds)  + (lmds-fmds)+1; % Data Set Sample Size @
    calds=zeros(nds,2); calds(1,1)=fyds; calds(1,2)=fmds; 
    yr=fyds; mt=fmds;
    i=2; for i=2:nds;
     mt=mt+1;
     if mt > 12; mt=1; yr=yr+1; end; 
     calds(i,1)=yr; calds(i,2)=mt;
     end; calds_init=calds;

    %Initial data
    fyreg=1959; fmreg=1; 
    cut=calds2n(calds,1959,1);
    y=y(cut:rows(y),:); m=m(cut:rows(m),:); calds=calds(cut:rows(calds),:); 
    gap=gap(cut:rows(gap),:); 
    y_sum=y_sum(cut:rows(y_sum),:);
    %Choose the variables and lags
    %im=1; %lags=max([lagsm;lagsy;h]); 
    R=calds2n(calds,yearR,monthR)-calds2n(calds,fyreg,fmreg)+1; 

    i_end=rows(m)-R; f_in_roll=[]; f_oos_roll=[]; series=[]; check=0; inflvar=[];
    bvec=[]; caldstoplotcoeff=[];
    for i=1:i_end-h;  
            t=R+i-1; t_init=R; 
            %Construct x and y
            x=[m,y];
            xf=[x]; 
              yhere=y(t+h-R+1+1:t,:)-y(t+h-R+1:t-1,:);  
              for s=1:h-1; yhere=yhere+y(t+h-R+2-s:t-s,:)-y(t+h-R+1-s:t-s-1,:);
              end; 
              x=x(t+h-R+1-h+1:t-h,:); 
            ty=length(yhere); ytrue=y(t+h)-y(t);
            biclags=bic(yhere,[y(t+h-R+1-h+1:t-h,:),m(t+h-R+1-h+1:t-h,:)],3,1); 
            x=[lagnmatrix0(m,biclags),lagnmatrix0(y,biclags)]; 
            xf=[ones(rows(x),1),lagnmatrix0(xf(:,1),biclags),lagnmatrix0(xf(:,2),biclags)];
            x=x(t+h-R+1-h+1:t-h,:); 
            yhere=yhere(biclags+1:rows(yhere)); 
            x=[ones(rows(x)-biclags,1), x(biclags+1:rows(x),:)]; 

            %Estimation
            results=ols(yhere,x); betahat=results.beta; 
            yhat=xf(t,:)*betahat;  
            f_in_roll = [f_in_roll ;  media((results.resid).^2) ];
            f_oos_roll= [f_oos_roll;  ( ytrue -  yhat)    .^2] ;
            series=[series;ytrue,yhat];

            %Estimation of the monetary policy rule
            regressors=[ones(rows(y),1),y,lagnmatrix0(ffr,2),gap];
            instruments=[ones(rows(y),1),lagnmatrix0(y,4),lagnmatrix0(gap,4),lagnmatrix0(ffr,4),lagnmatrix0(pi,4),lagnmatrix0(m2,4),lagnmatrix0(spread,4)];
            ffrhere=ffr(6+t-R:t,:); regressors=regressors(6+t-R:t,:); 
            instruments=instruments(6+t-R:t,:);
            rrstar=(sum(ffr)-sum(y))/rows(y); 
            b = gmmbeta(ffrhere,regressors,instruments,1); 
            b1 =[b(1),b(2)/(1-b(3)-b(4)),b(5)/(1-b(3)-b(4)),b(3)+b(4),rrstar,(rrstar -  b(1)/(1-b(3)-b(4))   )/((b(2)/(1-b(3)-b(4)))-1)];
            bvec=[bvec;b1]; 

            %Var of inflation
            inflvar=[inflvar;cov(yhere)];

            caldstoplotcoeff=[caldstoplotcoeff;calds(t,:)];

    end;%end loop for i

    lossdiff=f_oos_roll - f_in_roll;
    %disp('  '); disp(['Lags for unemployment = ',num2str(lagsm),' Lags for p =',num2str(lagsy)]);
    
        %Unconditional test
        n=rows(f_oos_roll);
        pai=rows(f_oos_roll)/R;
        if pai>1;
            [teststat ,pval ] = FF_test(lossdiff, 2/(3*pai)   );
        else;
            [teststat ,pval ] = FF_test(lossdiff, 1-(1/3)*(pai^2)   );
        end;
        disp(['Unconditional test p-value = ',num2str(pval)]);
        unconstat=teststat; pvalunc=pval; 

    
    %Conditional   test
    %constant only
        Zt=[ones(rows(lossdiff)-1,1) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_c,pval_FF_c]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),1,reshat,1);                
    disp(['Conditional test p-value (constant only) = ',num2str(pval_FF_c)]);
    %constant and lagged loss
        Zt=[ones(rows(lossdiff)-1,1),lossdiff(1:rows(lossdiff)-1,:) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_ar,pval_FF_ar]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0 1],reshat,1);                
    disp(['Conditional test p-value (constant and lagged loss) = ',num2str(pval_FF_ar)]);

    table=[table; h,R,pval,pval_FF_c,pval_FF_ar];

    %Using the monetary policy coefficients
    disp('Using the monetary policy coefficients');
    coefv=[]; sigmaev=h; 
    coefname=['const ';'beta  ';'gamma ';'rho   ';'rrstar';'pistar'];
    %Using coefficients one at a time -- marginal effect on mon pol coeff ONLY, NW HAC (h-1)
    pvstatv=[]; %%%pvstatv=[h,unconstat,FF_ar;h,pvalunc,pval_FF_ar]; 
        for j=2:4; 
            disp(coefname(j,:)); 
            Zt=[ones(rows(lossdiff),1),bvec(:,j)];
            coef=olsbeta(lossdiff,Zt),  reshat=olsres(lossdiff,Zt);
            [teststat,pval]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);
            pvstatv=[pvstatv,[coef(1,2);pval]];
        end;
    %Using inflation volatility
        Zt=[ones(rows(lossdiff)-h,1),inflvar(1:rows(lossdiff)-h,:)]; 
        coef=olsbeta(lossdiff(h+1:rows(lossdiff),:),Zt);
        reshat=olsres(lossdiff(h+1:rows(lossdiff),:),Zt);
        [teststat,pval]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);        
        pvstatv=[pvstatv,[coef(1,2);pval]];
    
   %All mon pol coeff simultaneously
        Zt=[ones(rows(lossdiff),1),bvec(:,2:4)];
        coef=olsbeta(lossdiff,Zt); reshat=olsres(lossdiff,Zt);       
        [teststat,pval_all]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[zeros(3,1),eye(3)],reshat,1);        
        pvstatv=[pvstatv,[teststat;pval_all]];
            
    Table3new=[Table3new; pvstatv(:,[1,2,3,5,4])];

    disp(['One-step ahead forecasts beginning at ',num2str(yearR),'  ',num2str(monthR)]);

    table, Table3new,
    diary off;

    end; %end loop yearmonthR
    
end;
diary off; 
delete TableBICnewer.out; diary TableBICnewer.out; disp('BIC RESULTS'); disp('Table2'); disp(table); disp('Table3'); disp(Table3new); diary off;