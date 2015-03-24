function omega = varPE(z,sigma,D,B,scheme,l,HAC_lag,pai); 
%This code provides the estimate omegahat(m,n) described in Corollary 8 of
%Giacomini and Rossi, December 2005. It is the estimated variance of the
%Wald test under conditional homoskedasticity. 
%The forecasting model is: SL(t+tau)=a+z(t)'*b+residual(t)
%where the SL are constructed based upon an estimated parameter
%betahat=inv(B)*H, where H is some sample average of h(t)'s

%INPUTS: 
%        z is the regressor WITHOUT THE CONSTANT
%        sigma is the variance of the out of sample losses
%        D is the sample mean of the derivative of the Surprise Loss relative to b
%        B see the formula above
%        l are the oos losses
%        'scheme' is either 1 (=fixed), 2(=rolling) or 3(=recursive)
%         pai is n/m
%         HAC_lag is the n of lags in the Newey and West HAC estimator
        
        l=l-mean(l,1); %demeaned losses
        if isempty(z)==1; ztilda=[]; else; 
        ztilda=z-ones(rows(z),1).*mean(z,1); end; %demeaned z
        Sll=nw(l,HAC_lag);
        Slzllzl=nw([l,mmult(l,ztilda)],HAC_lag); Slzl=Slzllzl(1,2:cols(Slzllzl));
        Szlzl=nw(mmult(l,ztilda),HAC_lag); 
        Szz=nw(ztilda,HAC_lag);
        if scheme==1; lambda=1; lambdahh=pai; lambdal=1+pai;
        elseif scheme==2; if pai<=1; lambda=1-pai/2; lambdahh=pai-(pai^2)/3; lambdal=1-(pai^2)/3; else;  lambda=1/(2*pai); lambdahh=1-1/(3*pai); lambdal=(2/3)/pai; end;
        elseif scheme==3; lambda=(log(1+pai))/pai; lambdahh=2*(1-(1/pai)*log(1+pai)); lambdal=1; end; 
        
        V11=sigma*lambdal;
        V12=lambda*Slzl;
        V22=Szlzl; 
        V=[V11,V12;V12',V22];
        
        if isempty(z)==1; MM=1; else; MM=[1,-mean(z,1)'*inv(Szz); zeros(cols(ztilda),1), inv(Szz)]; end;
        omega=MM*V*MM'; 
        