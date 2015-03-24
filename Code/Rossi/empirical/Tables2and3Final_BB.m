% This file replicates Tables 2 and 3 in Giacomini and Rossi,
% 'Detecting and Predicting Forecast Breakdown'
% LAGS = 1, 3 (Results for BIC are reported in Tables2and3BICFinal.m)
% REVISED DATA ONLY (Results for Real Time data are reported in Table2Final.m)
% Selected results are in Tables2and3.out
% This file also plots Figures 6 and 7
addpath('..\maincodes')
addpath('..\..\Functions')
clear; tableh=[]; ForecBreakd=[]; ForecBreakd1=[]; wanttoplot=1; Table3=[]; Table3new=[];
table=[]; ForecBreakd2=[]; 
hvec=12%[1,3,12]; %IF USING THIS PRG FOR TABLES
%hvec=[1:1:12]; %IF USING THIS PRG FOR PICTURES OF FORECAST LOSSES

for kk=1:cols(hvec); h=hvec(1,kk),

    yearmonthRv=[1979,1];%1985,1;1993,1];[1984,12; 1978,12; 1976,12; 1975,12; 1993,1]; 

    for yearmonthR=1:rows(yearmonthRv);
    yearR=yearmonthRv(yearmonthR,1); monthR=yearmonthRv(yearmonthR,2); 

    for ic=[1]; 
        for lagsm=1%1:2:3; 
            for lagsy=3%[1,3];
                lagsr=0; lagsp=0;

    cd data; cd qvmd; load ruc_latest.txt; load cpi_latest.txt; cd .. ; cd .. ;
    cd data; load fedfund5901.txt; load m2.txt; load commodp.txt; load tbill3m.txt; load tb10years.txt; cd .. ;

    y=cpi_latest; y_sum=[zeros(11,1);[100*log(y(13:rows(y),1)./y(1:rows(y)-12,1))]]; 
    y=1200*log(y(2:rows(y),:)./y(1:rows(y)-1,:)); y_init=y; 
    pi=commodp; pi=[0;1200*log(pi(2:rows(pi),:)./pi(1:rows(pi)-1,:))]; 
    m2=[0;log(m2(2:rows(m2),:)./m2(1:rows(m2)-1,:))];
    spread=tbill3m-tb10years; spread=[0;spread(2:rows(spread),:)];
    m=ruc_latest; m=m(2:rows(m),:); 
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
    im=1; lags=max([lagsm;lagsy;h]); 
    R=calds2n(calds,yearR,monthR)-calds2n(calds,fyreg,fmreg)+1; 

    i_end=rows(m)-R; f_in_roll=[]; f_oos_roll=[]; series=[]; check=0; inflvar=[];
    bvec=[]; caldstoplotcoeff=[]; f_all=[];
    f_in_btw=[];
    for i=1:i_end-h;  
            t=R+i-1; %t_init=R; 
            %Construct x and y
            if im==1; xm=lagnmatrix0(m,lagsm); x=[xm]; end;
            if lagsr>0; xr=lagnmatrix0(r,lagsr); x=[x,xr]; end;
            if ic==1; x=[x,ones(length(x),1)]; end;
            if lagsy>0; x=[x,lagnmatrix0(y,lagsy)]; end;
            if lagsp>0; x=[x,lagnmatrix0(p,lagsp)]; end; xf=x; 
              yhere=y(t+h-R+1+1:t,:)-y(t+h-R+1:t-1,:);  
              for s=1:h-1; yhere=yhere+y(t+h-R+2-s:t-s,:)-y(t+h-R+1-s:t-s-1,:);
              end; 
              x=x(t+h-R+1-h+1:t-h,:);  
            ty=length(yhere); ytrue=y(t+h)-y(t);

            %Estimation
            results=ols(yhere,x); betahat=results.beta; betahat_subsets=betahat(1:lagsm,:);
            yhat=xf(t,:)*betahat;  
            f_in_roll = [f_in_roll ;  media((results.resid).^2) ];
            f_oos_roll= [f_oos_roll;  ( ytrue -  yhat)    .^2] ;
            resul=(results.resid).^2; 
            if i==1; f_all=(results.resid).^2; end; 
            if h>1; if i<=h-1; f_all=[f_all;resul(rows(resul),1)]; end; end;
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

            %Calendar 
            caldstoplotcoeff=[caldstoplotcoeff;calds(t,:)];

    end;%end loop for i
    
    f_all=[f_all;f_oos_roll];
    
    lossdiff=f_oos_roll - f_in_roll;
    disp('  '); disp(['Lags for unemployment = ',num2str(lagsm),' Lags for p =',num2str(lagsy)]);

    %Unconditional test
    n=rows(f_oos_roll);
    pai=n/R;
    if pai>1;
        [teststat ,pval ] = FF_test(lossdiff, 2/(3*pai)   );
    else;
        [teststat ,pval ] = FF_test(lossdiff, 1-(1/3)*(pai^2)   );
    end;
    disp(['Unconditional test p-value = ',num2str(pval)]);
    unconstat=teststat; pvalunc=pval;
    [SL,SLbar,~,~,~,~,sigma2] = grtest(f_oos_roll,f_in_roll,y,xm,'rolling',round(n^(1/3)));

    %Conditional   tests  for the joint effect of the mon pol coeff and a constant (iid assn.)
    %constant only
    %[FF_c,pval_FF_c] = FF_test_cond(lossdiff(2:rows(lossdiff),:),); 
    %Zt=[ones(rows(lossdiff)-1,1) ]; 
    %coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt);
    %[FF_c,pval_FF_c]=FB_Waldtest(coef,Zt,f_all(2:rows(f_all),:),R-h-1,n-1,h,2,round(n^(1/3)),1);
        Zt=[ones(rows(lossdiff)-1,1) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_c,pval_FF_c]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),1,reshat,1);                
    disp(['Conditional test p-value (constant only) = ',num2str(pval_FF_c)]);

    %constant and lagged loss
    %Zt=[ones(rows(lossdiff)-1,1),lossdiff(1:rows(lossdiff)-1,:) ]; 
    %coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt);
    %[FF_ar,pval_FF_ar]=FB_Waldtest(coef,Zt,f_all(2:rows(f_all),:),R-h-1,n-1,h,2,round(n^(1/3)),1);
        Zt=[ones(rows(lossdiff)-1,1),lossdiff(1:rows(lossdiff)-1,:) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_ar,pval_FF_ar]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0 1],reshat,1);                
    disp(['Conditional test p-value (constant and lagged loss) = ',num2str(pval_FF_ar)]);
    
    table=[table;h,R,lagsm,lagsy,pval,pval_FF_ar];

    %Using the monetary policy coefficients
    disp('Using the monetary policy coefficients');
    coefv=[]; sigmaev=h; 
    coefname=['const ';'beta  ';'gamma ';'rho   ';'rrstar';'pistar'];
    
    %Using coefficients one at a time -- marginal effect on mon pol coeff ONLY, NW HAC (h-1)
    disp(['lagsm= ',num2str(lagsm),'; lagsy= ',num2str(lagsy)]); pvstatv=[h,lagsm,lagsy,unconstat,FF_ar;h,lagsm,lagsy,pvalunc,pval_FF_ar];
        for j=2:4; 
            disp(coefname(j,:)); 
            Zt=[ones(rows(lossdiff),1),bvec(:,j)];
            coef=olsbeta(lossdiff,Zt);  reshat=olsres(lossdiff,Zt);
            [teststat,pval]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);
            pvstatv=[pvstatv,[coef(1,2);pval]];    
            interv=FB_Waldtest_interval_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);
        figure; 
        plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,interv(:,1),'r','LineWidth',1.5);
        hold on; plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,interv(:,2),'g','LineWidth',1.5); 
        title([coefname(j,:),' at h ',num2str(h)]); hold off;
        end; 
 
    %Using inflation volatility
        Zt=[ones(rows(lossdiff)-h,1),inflvar(1:rows(lossdiff)-h,:)]; 
        coef=olsbeta(lossdiff(h+1:rows(lossdiff),:),Zt);
        reshat=olsres(lossdiff(h+1:rows(lossdiff),:),Zt);
        [teststat,pval]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);        
        pvstatv=[pvstatv,[coef(1,2);pval]];
        interv=FB_Waldtest_interval_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0,1],reshat,1);
        figure; 
        plot(caldstoplotcoeff(h+1:rows(caldstoplotcoeff),1)+caldstoplotcoeff(h+1:rows(caldstoplotcoeff),2)/12,interv(:,1),'r','LineWidth',1.5);
        hold on; plot(caldstoplotcoeff(h+1:rows(caldstoplotcoeff),1)+caldstoplotcoeff(h+1:rows(caldstoplotcoeff),2)/12,interv(:,2),'g','LineWidth',1.5); 
        title(['Inflation volatility',' at h ',num2str(h)]); hold off;

    %All mon pol coeff simultaneously
        Zt=[ones(rows(lossdiff),1),bvec(:,2:4)];
        coef=olsbeta(lossdiff,Zt); reshat=olsres(lossdiff,Zt);       
        [teststat,pval_all]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[zeros(3,1),eye(3)],reshat,1);        
        pvstatv=[pvstatv,[teststat;pval_all]];
        interv=FB_Waldtest_interval_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[zeros(3,1),eye(3)],reshat,1);
        [Wstat,pvalW,SLfitCI] = FB_Wald_CI(Zt,reshat,coef',n,R-h-1,sigma2,'rolling',round(n^(1/3)),0.05);
        figure; 
        plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,interv(:,1),'r','LineWidth',1.5);
        hold on; plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,interv(:,2),'g:','LineWidth',1.5);
        plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,interv(:,2).*0,'c','LineWidth',1.5);
        title('Fitted Surprise Losses and 95% confidence band'); 
        legend('Fitted Surprise Losses','95% Confidence Band'); hold off;
    Table3new=[Table3new;pvstatv(:,[1,2,3,6,7,8,10,9])]; 
    end; end; end;  

    disp(['One-step ahead forecasts beginning at ',num2str(yearR),'  ',num2str(monthR)]);

    table, Table3new, 
    diary off;

    end; %end loop yearmonthR

    %**********************************************************
    %PLOT FIGURE 6 -- rolling coefficient estimates over time
    %**********************************************************
    nn=rows(caldstoplotcoeff); nnv=[13:20:nn];
    plot(caldstoplotcoeff(nnv,1)+caldstoplotcoeff(nnv,2)/12,bvec(nnv,1),'k*','MarkerSize',10); hold on;
    plot(caldstoplotcoeff(nnv,1)+caldstoplotcoeff(nnv,2)/12,bvec(nnv,2),'ko','MarkerSize',10,'markerfacecolor','r'); 
    plot(caldstoplotcoeff(nnv,1)+caldstoplotcoeff(nnv,2)/12,bvec(nnv,3),'ks','MarkerSize',10,'markerfacecolor','m'); 
    plot(caldstoplotcoeff(nnv,1)+caldstoplotcoeff(nnv,2)/12,bvec(nnv,4),'kd','MarkerSize',10,'markerfacecolor','g'); 
    plot(caldstoplotcoeff(nnv,1)+caldstoplotcoeff(nnv,2)/12,bvec(nnv,5),'kv','MarkerSize',10,'markerfacecolor','b'); 
    legend('constant','\beta_t','\gamma_t','\rho_t','rr*');
    hold on; 
    plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,bvec(:,1),'k','LineWidth',1.5); hold on;
    plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,bvec(:,2),'r','LineWidth',1.5); 
    plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,bvec(:,3),'m','LineWidth',1.5); 
    plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,bvec(:,4),'g','LineWidth',1.5); 
    plot(caldstoplotcoeff(:,1)+caldstoplotcoeff(:,2)/12,bvec(:,5),'b','LineWidth',1.5); 
    axis tight; xlabel('Time'); 
    
    
    %*************************************************************
    %PLOTS FOR THE VARIOUS FORECAST SCENARIOS
    %*************************************************************
    %We first redo the estimation for the lags selected (here 3 and 3)
    %then we do the counterfactual analysis and plot the results
    
    %Using the monetary policy coefficients
    disp('Using the monetary policy coefficients');
    coefv=[]; sigmaev=h; 
    coefname=['const ';'beta  ';'gamma ';'rho   ';'rrstar';'pistar'];
    %Using coefficients one at a time -- marginal effect on mon pol coeff ONLY, NW HAC (h-1)
    disp(['lagsm= ',num2str(lagsm),'; lagsy= ',num2str(lagsy),' used for forecast scenarios']); pvstatv=[h;h];
        for j=2:4; %2
        disp(coefname(j,:)); 
        coef=[olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,j)]);olsstddevnw(lossdiff,[ones(rows(lossdiff),1),bvec(:,j)],h-1)'];
        teststat=olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,j)])./olsstddevnw(lossdiff,[ones(rows(lossdiff),1),bvec(:,j)],h-1)'; 
        pval= 1-chi2cdf(teststat.^2,1); disp([coef(1,:);pval]); pvstatv=[pvstatv,[coef(1,2);pval(1,2)]];
        res=olsres(lossdiff,[ones(rows(lossdiff),1),bvec(:,j)]); sigmae=sqrt(nw(res,h-1)); 
        end;
    %Using inflation volatility
        coef=([olsbeta(lossdiff,[ones(rows(lossdiff),1),inflvar]);olsstddevnw(lossdiff,[ones(rows(lossdiff),1),inflvar],h-1)']); 
        teststat=([olsbeta(lossdiff,[ones(rows(lossdiff),1),inflvar])./olsstddevnw(lossdiff,[ones(rows(lossdiff),1),inflvar],h-1)']); 
        pval= 1-chi2cdf(teststat.^2,1); 
        disp([coef(1,:);pval]); pvstatv=[pvstatv,[coef(1,2);pval(1,2)]];
       
    %All mon pol coeff simultaneously
        %olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)])
        coef=[olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)]);olsstddevnw(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)],h-1)'];
        teststat=olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)])./olsstddevnw(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)],h-1)';
        pval= 1-chi2cdf(teststat.^2,1); disp([coef(1,:);pval]);
        res=olsres(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)]); sigmae=sqrt(nw(res,h-1)); 
        
        %Various FORECAST SCENARIOS for BETA
        Lossforec=[]; coefv=[]; sigmaev=h; variousbetas=[0.83,2.15];
        for li=1:cols(variousbetas); betahere=variousbetas(1,li); 
        Lossforec=[1, bvec(rows(bvec),2),  betahere, bvec(rows(bvec),3:cols(bvec)-2)   ]*coef(1,:)'; %change rho
        coefv=[coefv,[olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)]);pval]];
        sigmaev=[sigmaev,Lossforec,sigmae];
        end; 
        
         %Various FORECAST SCENARIOS for RHO
        Lossforec=[]; coefv=[]; sigmaev1=h; variousrhos=[0.68,0.79,0.95,1,1.1];
        for ll=1:cols(variousrhos); rhohere=variousrhos(1,ll); 
        Lossforec=[1,  bvec(rows(bvec),1:cols(bvec)-3), rhohere   ]*coef(1,:)'; %change rho
        coefv=[coefv,[olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)]);pval]];
        sigmaev1=[sigmaev1,Lossforec,sigmae];
        end; 

        %Various FORECAST SCENARIOS for GAMMA
        Lossforec=[]; coefv=[]; sigmaev2=h; variousgammas=[0.27,0.93];
        for ll=1:cols(variousgammas); gammahere=variousgammas(1,ll); 
        Lossforec=[1,  bvec(rows(bvec),1:2), gammahere, bvec(rows(bvec),4)   ]*coef(1,:)'; %change rho
        coefv=[coefv,[olsbeta(lossdiff,[ones(rows(lossdiff),1),bvec(:,1:cols(bvec)-2)]);pval]];
        sigmaev2=[sigmaev2,Lossforec,sigmae];
        end; 
        
ForecBreakd=[ForecBreakd;sigmaev];
ForecBreakd1=[ForecBreakd1;sigmaev1];
ForecBreakd2=[ForecBreakd2;sigmaev2];
Table3=[Table3;pvstatv];
end;

for jj=1:round(cols(ForecBreakd)/2)-1;
figure; plot(ForecBreakd(:,1),ForecBreakd(:,jj*2),'LineWidth',1.8); hold on; 
plot(ForecBreakd(:,1),ForecBreakd(:,jj*2)+1.96*ForecBreakd(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd(:,1),ForecBreakd(:,jj*2)-1.96*ForecBreakd(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd(:,1),ForecBreakd(:,1)*0);
xlabel('\tau'); ylabel('SL_(_t_+_\tau_)'); hold off; title(['Estimated Surprise Losses, \beta=',num2str(variousbetas(1,jj))]);
end; 
    
for jj=1:round(cols(ForecBreakd1)/2)-1;
figure; plot(ForecBreakd1(:,1),ForecBreakd1(:,jj*2),'LineWidth',1.8); hold on; 
plot(ForecBreakd1(:,1),ForecBreakd1(:,jj*2)+1.96*ForecBreakd1(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd1(:,1),ForecBreakd1(:,jj*2)-1.96*ForecBreakd1(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd1(:,1),ForecBreakd1(:,1)*0)
xlabel('\tau'); ylabel('SL_(_t_+_\tau_)'); hold off; title(['Estimated Surprise Losses, \rho=',num2str(variousrhos(1,jj))]);
end; 

for jj=1:round(cols(ForecBreakd2)/2)-1;
figure; plot(ForecBreakd2(:,1),ForecBreakd2(:,jj*2),'LineWidth',1.8); hold on; 
plot(ForecBreakd2(:,1),ForecBreakd2(:,jj*2)+1.96*ForecBreakd2(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd2(:,1),ForecBreakd2(:,jj*2)-1.96*ForecBreakd2(:,jj*2+1),':','LineWidth',1.8);
plot(ForecBreakd2(:,1),ForecBreakd2(:,1)*0)
xlabel('\tau'); ylabel('SL_(_t_+_\tau_)'); hold off; title(['Estimated Surprise Losses, \gamma=',num2str(variousgammas(1,jj))]);
end; 


delete Tables2and3.out; diary Tables2and3.out; 
disp('Table 2 contains: const  R  q_u  q_pi  t_m,n,tau   W_m,n,tau_on_LaggedLosses');
disp('Table2'); disp(table); 
disp('Table3 contains: t_beta   t_gamma   t_rho   Wald_all_mon   t_infl_var');
disp('Table3'); disp(Table3new); diary off; 