# Sample length, training sample R and evaluation period P
R = 150
P = 150
T = R+P

# length of rolling window used to construct FLUC test
mu = 0.3
m = mu*P

# Number of simulations
Nsim = 5000

# Processes y and x
set.seed(1)
eps = matrix(rnorm((T*Nsim),mean = 0, sd = 1),T,Nsim)
nu = matrix(rnorm(T*Nsim,mean = 0, sd = 1),T,Nsim)

x = matrix(0,T,Nsim)
x[1,] = nu[1,]
for (t in 2:T){
  x[t,] = 0.5*x[(t-1),]+nu[t,]
}
beta = 0.1/sqrt(R)
y = beta*x+eps

# Estimate models and obtain forecasts and forecast errors
# Forecasts f1 and f2
f1 = matrix(0,P,Nsim)
f2 = matrix(0,P,Nsim)
# Squared forecast errors l1 and l2
l1 = matrix(0,P,Nsim)
l2 = matrix(0,P,Nsim)
for (n in 1:Nsim){
  for (j in R:(T-1)){
    xest = x[((j-R+1):j),n]
    yest = y[((j-R+1):j),n]
    b = solve(t(xest)%*%xest)%*%t(xest)%*%yest
    f2[(j-R+1),n] = b%*%x[(j+1),n]
    l1[(j-R+1),n] = y[(j+1),n]^2
    l2[(j-R+1),n] = (y[(j+1),n]-f2[(j-R+1),n])^2
  }
}

# Loss differential
ldiff = l1-l2

# FLUC Test
resultsfluc = vector('list',Nsim)
testsizefluc = 0
for (n in 1:Nsim){
  resultsfluc[[n]] = fluc(ldiff[,n,drop=F],mu)
  testsizefluc = testsizefluc+resultsfluc[[n]]$rej
}
testsizefluc = testsizefluc/Nsim

# One-time reversal test
resultsotrev = vector('list',Nsim)
testsizeotrev = 0
for (n in 1:Nsim){
  resultsotrev[[n]] = otreverse(ldiff[,n,drop=F])
  testsizeotrev = testsizeotrev+resultsotrev[[n]]$rej
}
testsizeotrev = testsizeotrev/Nsim