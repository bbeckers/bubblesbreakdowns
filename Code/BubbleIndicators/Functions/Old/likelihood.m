function lik = likelihood(y,x,beta,sigma2)

T = size(y,1);

lik = -(T/2)*log(2*pi*sigma2)-sum((y-x'*beta)^2)/(2*sigma2);