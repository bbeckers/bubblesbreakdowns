function [beta,F] = mz_test(Error,Yfc,h)

beta = zeros(2*length(h),2);
F = zeros(length(h),2);
hind = 1;
for j=h'
    y = -(Error(hind,~isnan(Error(hind,:))))';
    T = size(y,1);
    X = [ones(T,1) Yfc(hind,1:T)'];
    betaj = (X'*X)^(-1)*X'*y;
    ehat = y-X*betaj;
    sigma2 = (ehat'*ehat)/(T-2);
    if j>1
        sigma2j = zeros(j-1,1);
        weights = zeros(j-1,1);
        for i=1:j-1
            sampleCov = cov(ehat(1+i:T),ehat(1:T-i));
            sigma2j(i) = sampleCov(2);
            weights(i) = 1-i/(j-1+1);
        end
        sigma2 = sigma2+2*sum(weights.*sigma2j);
    end
    sigma2_beta = sigma2*(X'*X)^(-1);
    t_beta = (betaj-zeros(2,1))./sqrt(diag(sigma2_beta));
    p_beta = 2*(1-tcdf(abs(t_beta),T-2));
    F(hind,1) = (betaj-zeros(2,1))'*sigma2_beta^(-1)*(betaj-zeros(2,1))/2;
    F(hind,2) = 1-fcdf(F(2*hind-1),2,T-2);
    beta(2*hind-1:2*hind,:) = [betaj'; p_beta'];
    hind = hind+1;
end