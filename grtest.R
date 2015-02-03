grtest=function(LO,LI,Y,X,tau,s='fixed'){        
        ## Description ##
        # ----------------------------------------------------------------------- #
        # This function conducts the forecast breakdown test by Giacomini & Rossi
        # (2009, The Review of Economic Studies).
        # ----------------------------------------------------------------------- #
        # INPUT to the function:
        # LO:  [n] Matrix of out-of-sample losses
        # LI:  [n] Matrix of average in-sample losses for each estimation point t
        # Y:   [T] Matrix of target variable
        # X:   [(T)x(k)] Matrix of predictors
        # m:   number of obs. for estimation ()
        # n:   number of in- and out-of-sample losses considered (fixed window)
        # T:   m+n
        # tau: Scalar denoting forecast horizon
        # s:   String sececting the forecasting scheme
        #      'fixed' (default), 'rolling', 'recursive'
        # ----------------------------------------------------------------------- #
        # OUTPUT of the function:
        # SL:    [n] Matrix of surprise losses
        # SLbar: Scalar average surprise loss
        # tc:    Overfitting-corrected test-statistic
        # pc:    Overfitting-corrected p-value
        # t:     Uncorrected test statistic
        # p:     Uncorrected p-value
        # ----------------------------------------------------------------------- #
        
        ## Begin function
        T = nrow(X)
        k=ncol(X)
        n = nrow(LO)
        m = T-n-tau+2
        bw = n^(1/3)
        # bw = 0
        
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
        else{SLL = 0}
        for (j in 1:bw){
                Gamma = t(LOtilde[(1+j):n])%*%LOtilde[1:(n-j),1]/nrow(LOtilde)
                SLL = SLL+2*(1-j/(bw+1))*Gamma
        }
        SLL = SLL+t(LOtilde)%*%LOtilde/nrow(LOtilde)
        
        # Variance estimator
        sigma = sqrt(lambda*SLL)
        
        ## Test statistic
        t = sqrt(n)*SLbar/sigma
        p =2*(1-pnorm(abs(t)))
        
        
        ## Overfitting correction
        # Asymptotic variance covariance matrix of coefficients
        
        XX = t(X[1:(T-tau),])%*%X[1:(T-tau),]
        
        beta = solve(XX)%*%t(X[1:(T-tau),])%*%Y[(1+tau):T,1]
        u = Y[(1+tau):T,1]-X[1:(T-tau),]*beta
        Vbeta = t(u)%*%u/(T-tau-k)%*%solve(XX)
        # Correction parameter
        
        c = 2*gamma*tr((T-tau)^(-1)*XX%*%Vbeta)
        tc = (sqrt(n)*SLbar-c)/sigma
        pc = 2*(1-pnorm(abs(t)))
        results=list(c(SL,SLbar,tc,pc,t,p))
        return(results)
}# function brackets