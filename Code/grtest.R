grtest=function(LO,LI,Y,X,scheme,bw){        
        ## Description
        # ----------------------------------------------------------------------- #
        # This function conducts the forecast breakdown test by Giacomini & Rossi
        # (2009, The Review of Economic Studies).
        # ----------------------------------------------------------------------- #
        # INPUT to the function:
        # LO:  [n] Vector of out-of-sample losses
        # LI:  [n] Vector of average in-sample losses for each estimation point t
        # Y:   [T] Vector of target variable
        # X:   [(T)x(k)] Matrix of predictors
        # tau: Scalar denoting forecast horizon
        # s:   String sececting the forecasting scheme
        #      'fixed' (default), 'rolling', 'recursive'
        # ----------------------------------------------------------------------- #
        # OUTPUT of the function:
        # SL:     [n] Vector of surprise losses
        # SLbar:  Scalar average surprise loss
        # tc:     Overfitting-corrected test-statistic
        # pc:     Overfitting-corrected p-value
        # t:      Uncorrected test statistic
        # p:    Uncorrected p-value
        # sigma2: Variance of out-of-sample surprise losses
        # ----------------------------------------------------------------------- #
        
        ## Begin function
        T = nrow(X)
        k = ncol(X)
        n = nrow(LO)
        m = T-n
        
        # Adjustment paramater for surprise loss variance estimator
        lambda=1
        if (s=='fixed'){lambda = 1+n/m}
        if (s=='rolling'& n<m){lambda = 1-(1/3)*(n/m)^2}
        if (s=='rolling'&n>=m){lambda = (2/3)*(m/n)}
        
        # Adjustment paramater in overfitting correction
        if (s=='recursive'){gamma = (1/sqrt(n))*log(1+n/m)}
        else{gamma = sqrt(n)/m}
        
        ## Average surprise losses
        # Matrix of surprise losses
        SL = LO-LI
        
        # Out-of-sample mean of surprise losses
        SLbar = mean(SL)
        
        ## Asymptotic variance
        # Demeaned surprise losses
        LOtilde = LO-mean(LO)
        
        # HAC variance estimator of demeaned surprise losses
        if (bw==0){SLL = cov(LO)}
        else{SLL = 0
          for (j in 1:bw){
                Gamma = t(LOtilde[(1+j):n])%*%LOtilde[1:(n-j),1]/n
                SLL = SLL+2*(1-j/(bw+1))*Gamma
          }
        SLL = SLL+t(LOtilde)%*%LOtilde/n
        }
        
        # Variance estimator
        sigma2 = lambda*SLL
        
        ## Test statistic
        t = sqrt(n)*SLbar/sqrt(sigma2)
        p = 1-pnorm(abs(t))
        
        
        ## Overfitting correction
        # Asymptotic variance covariance matrix of coefficients
        
        XX = t(X)%*%X
        
        beta = solve(XX)%*%t(X)%*%Y
        u = Y-X%*%beta
        Vbeta = (t(u)%*%u)[1,1]/(T-k)*solve(XX)
        # Correction parameter
        
        c = 2*gamma*sum(diag(XX%*%Vbeta))
        tc = t-c/sqrt(sigma2)
        pc = 1-pnorm(abs(t))
        results = list(SL,SLbar,tc,pc,t,p,sigma2)
        return(results)
}# function brackets