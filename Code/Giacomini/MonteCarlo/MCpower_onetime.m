% MC experiment 

% size: dgp is y_t=phi*x_t+eps_t, eps_t iid N(0,sig2), where x_t is iid N(0,1)
% (drawn once and kept fixed in all experiments),
% competing models are y_t=eps_t and y_t=phi*x_t+eps_t,
% for a pair of in-sample and out-of-sample sizes (R,n), models perform
% equally well at all points in time when
% phi=sig*sqrt(sum(x_t+1^2/Sxx-n)/sum(x_t+1^2))

function [nrejfluc,nrejDM,nrejopt]=MCpower_onetime(R,n,m,del)



nrep = 1000;
rej1 = zeros(nrep,1);
rej2 = zeros(nrep,1);
rej3 = zeros(nrep,1);
sig = 1;
T = R + n;



bb1 = sig/sqrt(R)-del;


bb2 = sig/sqrt(R)+del;

btot = [bb1*ones(R+(2*n)/3,1);bb2*ones(n/3,1)];
for rep = 1:nrep

x=randn(T,1);
for tt=2:T
    x(tt)=.5*x(tt-1)+randn;
end
    y = btot.*x + sig.*randn(T,1);
% num=0;
% den=0;
%      for ii=R:T-1
%     xest = x(ii-R+1:ii);
% 
%     sxx = sum(xest.^2);
%     num=num+(1/sxx)*x(ii+1)^2;
%     den = den+x(ii+1)^2;
% end
% bb=sig*sqrt(num/den);

    % construct oos forecast loss sequences from two models

    l1 = zeros(n,1);
    l2 = zeros(n,1);
    jj = 1;
    for ii = R:T-1
        f1 = 0;
        xest = x(ii-R+1:ii);
        yest = y(ii-R+1:ii);
        best=xest\yest;
%    best
        f2 = best*x(ii+1);
        l1(jj) = (y(ii+1)-f1)^2;
        l2(jj) = (y(ii+1)-f2)^2;
        jj = jj +1;
    end
% [    mean(l1) sig^2/R ]
% [    mean(l2) bb1^2 ]

    ldiff = l1-l2;

    % do the test
    fluc = zeros(n-m,1);
    jj = 1;
    % variance for fluctuation statistic computed over oos portion
    sigmahat2 = var(ldiff);
    for ii=m/2+1:n-m/2
        ldiffest = ldiff(ii-m/2:ii+m/2);
        fluc(jj) = (1/sqrt(sigmahat2*m))*sum(ldiffest);
        jj = jj + 1;
    end

    mu = m/n;
    if mu ==.1
        crit = 3.39;
    elseif mu == .2
        crit = 3.17;
    elseif mu == .3
        crit = 3.02;
    elseif mu == .4
        crit = 2.89;
    elseif mu == .5
        crit = 2.77;
    elseif mu == .6
        crit = 2.63;
    elseif mu == .7
        crit = 2.56;
    elseif mu == .8
        crit = 2.43;
    elseif mu == .9
        crit = 2.24;
    else
        'not integer'
    end
    if max(abs(fluc))>crit
        rej1(rep) =  1;
    end
    DM = sqrt(n)*mean(ldiff)/std(ldiff);
    if abs(DM)>1.96
        rej2(rep) = 1;
    end
    low = round(0.15*n);
    up = round(0.85*n);
    sig2 = var(ldiff);
    sizew = up-low;
    lm2=zeros(sizew,1);
    lm1=((sum(ldiff))^2)/(sig2*n);
    jj = 1;
    for ii=low:up
        a = 1/(sig2*ii);
        b = 1/(1-ii/n);
        lm2(jj) = a*b*(sum(ldiff(1:ii))-(ii/n)*sum(ldiff))^2;
        jj = jj + 1;
    end

    tstat = max(lm1+lm2);
    if tstat>9.82
        rej3(rep) =  1;
    end
end
nrejfluc = mean(rej1);
nrejDM = mean(rej2);
nrejopt = mean(rej3);



