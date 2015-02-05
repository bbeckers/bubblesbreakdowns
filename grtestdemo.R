# demo of giacomini rossi 2009 code
T=400
Tstar=300


Nsim = 1

# Parameters
alphax = 3
betax = 0.5
alphay1 = 1.3
betay1 = 0.2
alphay2 = 2.8
betay2 = 0.7

# Initial training sample size
m0 = 60
# Evaluation window
n = 40
# Forecast horizon
tau = 1
# Bandwidth for HAC estimator
bw = 0

# List (formerly cell) of regressors of insample regression

XY=list()

# Estimation scheme
s = 'recursive'
n.bt=T-m0+1-tau-n# number of breakdown tests
# Output
SL = array(NA,dim=c(Nsim,n,n.bt))
SLbar = matrix(NA,n.bt,Nsim) # mean surprise losses
tc = matrix(NA,n.bt,Nsim) # t stat with overfitting correction
pc = matrix(NA,n.bt,Nsim) # p value with overf. corr.
t = matrix(NA,n.bt,Nsim) # t stat simple
p = matrix(NA,n.bt,Nsim) # p stat simple

tgr = matrix(NA,n.bt,Nsim) # giacomini rossi from original code
pgr = matrix(NA,n.bt,Nsim) # p value of giac. ross. from orig. code
tgrc = matrix(NA,n.bt,Nsim) # giacomini rossi from original code (overfitt. corr.)
pgrc = matrix(NA,n.bt,Nsim)# p value of giac. ross. from orig. code (overfitt. corr.)

# Simulate time series x (predictor) and y (target)

for (nsim in 1:Nsim){# nsim=1
        set.seed(nsim)
        epsx = matrix(rnorm(T+1,0,1),ncol=1) # white noise
        epsy = matrix(rnorm(T,0,1),ncol=1) # white noise
        results=vector('list',n.bt)
        # Time series with break in y
        x = matrix(NA,T+1,1)
        
        y = matrix(NA,T,1)
        x[1,1] = alphax
        for (j in 2:(T+1)){
                x[j,1] = alphax+betax*x[j-1,1]+epsx[j,1]
        }
        
        for (j in 1:Tstar){
                y[j,1] = alphay1+betay1*x[j,1]+epsy[j,1]
        }
        
        for (j in (Tstar+1):T){
                y[j,1] = alphay2+betay2*x[j,1]+epsy[j,1]
        }
        rm(j)
        
        # Inputs to the function
        LO = matrix(NA,n,n.bt)
        LI = matrix(NA,n,n.bt)
        X=matrix(NA,T,n.bt)
        Y=matrix(NA,T,n.bt)
                
        # Begin forecast
        for (m in m0:(T-tau-n)){# m=m0 ; computing n-prediction errors.
                for (l in (1:n)){# l=1
                        Yest = y[1:(m+l-1),1]
                        ones = matrix(1,nrow=(m+l-1),ncol=1)
                        Xest = cbind(ones,x[1:(m+l-1),1])
                        beta = solve((t(Xest)%*%Xest))%*%t(Xest)%*%Yest
                        Yfit = Xest%*%beta
                        u = Yfit-Yest
                        LI[l,m-m0+1] = mean(u^2)
                        yfc = matrix(c(1,x[m+l-1]),nrow=1,ncol=2)%*%beta
                        e = yfc-y[m+l]
                        LO[l,m-m0+1] = e^2
                }
                Y[1:(m+l-1),m-m0+1] = Yest
                X[1:(m+l-1),m-m0+1] = Xest[,2]
                #                 result = linear_FF(Yest,Xest,m,tau,'sequ')
                #                 resultc = linear_FF_c(Yest,Xest,m,tau,'sequ')
                #                 tgr[m-m0+1,nsim] = result(1)        
                #                 pgr[m-m0+1,nsim] = result(2)
                #                 tgrc[m-m0+1,nsim] = resultc(1)        
                #                 pgrc[m-m0+1,nsim] = resultc(2)
        }
        
        # Evaluate forecast performance
        for (j in 1:(T-tau-n-m0+1)){
                nna=is.na(Y[,j])==F
                results[j] = grtest(LO[,j,drop=F],LI[,j,drop=F],Y[nna,j,drop=F],cbind(1,X[nna,j,drop=F]),s,bw)
        }
        nsim
}