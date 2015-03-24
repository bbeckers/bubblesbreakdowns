function omega = varPE(z,sigma,D,B,h,scheme,l,HAC_lag,pai); 
%The forecasting model is: SL(t+tau)=a+z(t)'*b+residual(t)
%where the SL are constructed based upon an estimated parameter
%betahat=inv(B)*H
%where H is some sample average of h(t)'s

%INPUTS: 
%        z is the regressor WITHOUT THE CONSTANT
%        sigma is the variance of the out of sample losses
%        D is the sample mean of the derivative of the Surprise Loss relative to b
%        B, h see the formula above
%        l are the oos losses
%        'scheme' is either 1 (=fixed), 2(=rolling) or 3(=recursive)
%         pai is n/m
%         HAC_lag is the n of lags in the Newey and West HAC estimator
        
        %H=mean(h,1); 
        l=l-mean(l,1); %demeaned losses
        if isempty(z)==1; ztilda=[]; else; 
        ztilda=z-ones(rows(z),1).*mean(z,1); end; %demeaned z
        %S=nw([l,h],HAC_lag); Sll=S(1,1); Slh=S(1,2:cols(S)); Shh=S(2:cols(S),2:cols(S));
        Sll=nw(l,HAC_lag);
        Slzllzl=nw([l,mmult(l,ztilda)],HAC_lag); Slzl=Slzllzl(1,2:cols(Slzllzl));
        %Shzlhzl=nw([h,mmult(l,ztilda)],HAC_lag); Shzl=Shzlhzl(1,2:cols(Slzllzl));
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
        