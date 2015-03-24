function [teststat,pvalue] = FB_Waldtest(deltahat,Z,losses_all,m,n,tau,scheme,hac_lag,R)

% Performs the Wald test discussed in Proposition 8 of Giacomini and Rossi (2005)
% for a quadratic loss.
%
%
% OUTPUT: pvalue, the p-value of the Wald test
%         omegahat, the estimator of the variance of deltahat
% INPUTS: deltahat, the OLS coefficient estimate in the regression of the surprise losses on explanatory
%         variables Z
%         Z, the nxr matrix of explanatory variables (typically including a
%         column of 1)
%         losses_all, a Tx1 vector of stacked in-sample and out-of-sample losses
%                    (see Algorithm 2 in the paper for explanation), where
%                    T = m + n + tau - 1
%         m, the in-sample size
%         n, the out-of-sample size
%         tau, the forecast horizon
%         scheme, the forecasting scheme. scheme = 1 for fixed; scheme = 2
%                 for rolling; scheme = 3 for recursive
%         hac_lag, the truncation lag for the Newey-West estimator (a
%                  rule-of-thumb choice is n^(1/3)) of the asymptotic
%                  variance
%         R, the matrix of restrictions (dimension #restr*r)


T = m + n + tau - 1;
r = size(Z,2);
lossdem = losses_all-mean(losses_all); 

% -------------------------------- FIXED SCHEME --------------------

if scheme == 1
    % weighted losses
    wl = zeros(T,r);
    for ii=1:m
        wl(ii,:) =  -(sum(Z)/m)*lossdem(ii);
    end
    for ii=m+tau:T
        wl(ii,:) = Z(ii-m-tau+1,:)*lossdem(ii);
    end
    % Compute the Newey-West estimator of the variance of the weighted losses
    wl = wl - ones(T,1)*mean(wl);
    samplevar = wl'*wl/T; % sample variance
    gamma = -999*ones(hac_lag,r);
    if hac_lag > 0
        % sample autocovariances
        for ii = 1:hac_lag
            lag = [zeros(ii,r);wl(1:T-ii,:)];
            gamma = (wl'*lag +lag'*wl)/T;
            weights = 1 - (ii/(hac_lag+1));
            samplevar = samplevar + weights*gamma;
        end
    end
    omegahat =  T*inv(Z'*Z)*samplevar*inv(Z'*Z);


    % -------------------------------- ROLLING SCHEME --------------------

elseif scheme == 2
    if n < m
        % weighted losses
        wl = zeros(T,r);
        for ii=1:n
            wl(ii,:) =  -(sum(Z(1:ii,:))/m)*lossdem(ii);
        end
        for ii=n+1:m
            wl(ii,:) =  -(sum(Z)/m)*lossdem(ii);
        end
        if tau>1
            for ii=m+1:m+tau-1
                wl(ii,:) =  -(sum(Z(ii-m+1:end,:))/m)*lossdem(ii);
            end
        end
        for ii=m+tau:T-tau
            wl(ii,:) =  (Z(ii-m-tau+1,:)-(sum(Z(ii-m+1:end,:))/m))*lossdem(ii);
        end
        for ii=T-tau+1:T
            wl(ii,:) =  Z(ii-m-tau+1,:)*lossdem(ii);
        end
        % Compute the Newey-West estimator of the variance of the weighted losses
        wl = wl - ones(T,1)*mean(wl);
        samplevar = wl'*wl/T; % sample variance
        gamma = -999*ones(hac_lag,r);
        if hac_lag > 0
            % sample autocovariances            
            for ii = 1:hac_lag
                lag = [zeros(ii,r);wl(1:n-ii,:)];  
                gamma = (wl'*lag +lag'*wl)/n;
                weights = 1 - (ii/(hac_lag+1));
                samplevar = samplevar + weights*gamma;
            end
        end
        omegahat =  T*inv(Z'*Z)*samplevar*inv(Z'*Z);
    else
        % weighted losses
        wl = zeros(T,r);
        for ii=1:m
            wl(ii,:) =  -(sum(Z(1:ii,:))/m)*lossdem(ii);
        end
        if tau>1
            for ii=m+1:m+tau-1
                wl(ii,:) =  -(sum(Z(ii-m+1:ii,:))/m)*lossdem(ii);
            end
        end
        for ii=m+tau:n
            wl(ii,:) =  (Z(ii-m-tau+1,:)-(sum(Z(ii-m+1:ii,:))/m))*lossdem(ii);
        end
        for ii=n+1:T-tau
            wl(ii,:) =  (Z(ii-m-tau+1,:)-(sum(Z(ii-m+1:end,:))/m))*lossdem(ii);
        end
        for ii=T-tau+1:T
            wl(ii,:) =  Z(ii-m-tau+1,:)*lossdem(ii);
        end
        % Compute the Newey-West estimator of the variance of the weighted losses
        wl = wl - ones(T,1)*mean(wl); 
        samplevar = wl'*wl/T; % sample variance
        samplevar=cov(wl);
        gamma = -999*ones(hac_lag,r); 
        if hac_lag > 0
            % sample autocovariances
            for ii = 1:hac_lag
                lag = [zeros(ii,r);wl(1:T-ii,:)]; 
                gamma = (wl'*lag +lag'*wl)/T;
                weights = 1 - (ii/(hac_lag+1));
                samplevar = samplevar + weights*gamma;
            end
        end
        omegahat =  T*inv(Z'*Z)*samplevar*inv(Z'*Z);
    end


    % -------------------------------- RECURSIVE SCHEME --------------------

elseif scheme == 3
    % weighted losses
    wl = zeros(T,r);
    am = zeros(r,n);
    for ii = 1:n
        for jj=ii-1:n-1
            am(:,ii) = am(:,ii)+Z(jj+1,:)'/(m+jj);
        end
    end
    for ii=1:m
        wl(ii,:) =  -am(:,1)'*lossdem(ii);
    end
    if tau>1
        for ii=m+1:m+tau-1
            wl(ii,:) =  -am(:,ii-m+1)'*lossdem(ii);
        end
    end
    for ii=m+tau:T-tau
        wl(ii,:) =  (Z(ii-m-tau+1,:)-am(:,ii-m+1)')*lossdem(ii);
    end
    for ii=T-tau+1:T
        wl(ii,:) =  Z(ii-m-tau+1,:)*lossdem(ii);
    end
    % Compute the Newey-West estimator of the variance of the weighted losses
    wl = wl - ones(T,1)*mean(wl);
    samplevar = wl'*wl/T; % sample variance
    gamma = -999*ones(hac_lag,r);
    if hac_lag > 0
        % sample autocovariances
        for ii = 1:hac_lag
            lag = [zeros(ii,r);wl(1:T-ii,:)];
            gamma = (wl'*lag +lag'*wl)/T;
            weights = 1 - (ii/(hac_lag+1));
            samplevar = samplevar + weights*gamma;
        end
    end
    omegahat =  T*inv(Z'*Z)*samplevar*inv(Z'*Z);
end

%==================================================================
% Wald test

teststat = (R*deltahat')'*inv(R*omegahat*R')*(R*deltahat');
pvalue = 1 - cdf('chi2',teststat,rows(R));
