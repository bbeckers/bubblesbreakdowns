%THIS IS THE RIGHT FILE TO DO EMPIRICAL FOR REALTIME DATA for BIC results
clear;  table=[];
cd data; load ruc7801.txt; cd .. ;
cd newdata; load cpi_sa_nan_jan78; cd .. ; 
y=cpi_sa_nan_jan78;  
y_sum=[zeros(12,rows(cpi_sa_nan_jan78));zeros(rows(cpi_sa_nan_jan78)-12,12),[100*log(y(13:rows(y),13:cols(y))./y(1:rows(y)-12,13:cols(y)))]]; 
%y_sum=[zeros(12,rows(cpi_sa_nan_jan78));zeros(rows(cpi_sa_nan_jan78)-12,12),[(y(13:rows(y),13:cols(y)))]]; TO DOUBLE CHECK TIMING ONLY

y=[zeros(1,rows(cpi_sa_nan_jan78));[zeros(rows(cpi_sa_nan_jan78)-1,1),1200*log(y(2:rows(y),2:cols(y))./y(1:rows(y)-1,2:cols(y)))]];
%y=[zeros(1,rows(cpi_sa_nan_jan78));[zeros(rows(cpi_sa_nan_jan78)-1,1),(cpi_sa_nan_jan78(2:rows(cpi_sa_nan_jan78),2:cols(cpi_sa_nan_jan78)))]]; %TO DOUBLE CHECK TIMING ONLY
y=y(4:rows(y),4:cols(y)); y_sum=y_sum(4:rows(y_sum),4:cols(y_sum)); 
m=ruc7801;  m=m(4:rows(m),2:cols(m)); 

%Data set information
fyds=1978; fmds=4;                                  % First Year (fyds) and quarter (fqds) of Data Set and vintage (fyvt)%
lyds=2001; lmds=12;                                 % Last Year and quarter of Data Set %
calds=calendarm(fyds,fmds,lyds,lmds);               % Set calendar for monthly
plot(calds(13:rows(calds),1)+calds(13:rows(calds),2)./12,diag(y_sum(13:rows(y_sum),13:rows(y_sum))))
title('Inflation with real time data'); axis tight;

% Select only quarterly variables
q_index=1:3:rows(calds);                                
m=m(q_index,:); y=y(q_index,q_index); calds=calds(q_index,:); y_sum=y_sum(q_index,q_index); 
%disp('The only possible dates to start the regression are: ');
%disp('Recall that this is the day of the data -- vintage is that plus one month');
%disp(calds); 
yearR=1993; monthR=1; lyearR=2001; lmonthR=10;      % First and Last estimation period
R=calds2n(calds,yearR,monthR)-calds2n(calds,fyds,fmds)+1; 

%Choose the variables and lags
for h=[1,4];   

    %for lagsm=lagsmv; for lagsy=lagsyv; lags=max([lagsm;lagsy;h]);

        f_in_roll=[]; f_oos_roll=[]; ytruev=[];
        for i=1:(rows(m)-h-R+1);   

                t=R+i-1; t_init=R;

                %Construct x and y
                x=[m(:,t),y(:,t)];
                xf=[x]; lags=h;
                yhere=y_sum(lags+1+1:t,t);
                %x=x(lags+1-h+1+1:t-h+1,:);  %t-h +ONE b/o x was already lagged once
                ty=length(yhere); ytrue=y_sum(t+h,t+h);   
                biclags=bic(yhere,[y(lags+1-h+1+1:t-h+1,t),m(lags+1-h+1+1:t-h+1,t)],3,1); %t+h-R+1-h+1:t-h,t
                x=[lagnmatrix0(m(:,t),biclags),lagnmatrix0(y(:,t),biclags)]; 
                xf=[ones(rows(x),1),lagnmatrix0(xf(:,1),biclags),lagnmatrix0(xf(:,2),biclags)];
                x=x(lags+1-h+1+1:t-h+1,:); 
                ytruev=[ytruev;ytrue]; 
                x=[ones(rows(x),1), x]; 

                %Estimation
                results=ols(yhere,x); betahat=results.beta; %betahat_subsets=betahat(1:lagsm,:);
                yhat=xf(t,:)*betahat;  
                f_in_roll = [f_in_roll ;  media((results.resid).^2) ];
                f_oos_roll= [f_oos_roll;  ( ytrue -  yhat)    .^2] ;

        end;%end loop for i

        lossdiff=f_in_roll - f_oos_roll;
        disp('  ');
        %disp(['Lags for unemployment = ',num2str(lagsm),' Lags for p =',num2str(lagsy)]);
        
        %Unconditional test
        n=rows(f_oos_roll);
        pai=rows(f_oos_roll)/R;
        if pai>1;
            [teststat ,pval ] = FF_test(lossdiff, 2/(3*pai)   );
        else;
            [teststat ,pval ] = FF_test(lossdiff, 1-(1/3)*(pai^2)   );
        end;
        disp(['Unconditional test p-value = ',num2str(pval)]);
        
        %Conditional   test
        %constant only
        Zt=[ones(rows(lossdiff)-1,1) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_c,pval_FF_c]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),1,reshat,1);                
        %constant and lagged loss
        Zt=[ones(rows(lossdiff)-1,1),lossdiff(1:rows(lossdiff)-1,:) ];
        coef=olsbeta(lossdiff(2:rows(lossdiff),:),Zt); reshat=olsres(lossdiff(2:rows(lossdiff),:),Zt);
        [FF_ar,pval_FF_ar]=FB_Waldtest_covstat(coef,Zt,R-h-1,n,h,2,round(n^(1/3)),[0 1],reshat,1);                
        disp(['Conditional test p-value (constant and lagged loss) = ',num2str(pval_FF_ar)]);

        table=[table;h,pval,pval_FF_ar];

end;  %closes the loop for h

delete Table2RealTimeBICFinal.out; diary Table2RealTimeBICFinal.out;
disp(' ');
disp('h, lagsm, lagsy, pval, pval_FF_ar');
disp(table);
diary off; 