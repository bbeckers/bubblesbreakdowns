function samplevar = NeweyWest(wl,hac_lag); 

%% Compute the Newey-West estimator
    T=rows(wl); 
    wl = wl - ones(T,1)*mean(wl);
    samplevar = wl'*wl/T; % sample variance
    %gamma = -999*ones(hac_lag,r);
    if hac_lag > 0
        % sample autocovariances
        for ii = 1:hac_lag
            lag = [zeros(ii,cols(wl));wl(1:T-ii,:)];
            gamma = (wl'*lag +lag'*wl)/T;
            weights = 1 - (ii/(hac_lag+1));
            samplevar = samplevar + weights*gamma;
        end
    end
