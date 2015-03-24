function [hADF_notrend,hADF_trend,pValADF_notrend,pValADF_trend,...
    hKPSS_notrend,hKPSS_trend,pValKPSS_notrend,pValKPSS_trend,...
    pstarnotrend,pstartrend,trend,trendpVal,hADF,pValADF,hKPSS,pValKPSS] = stationarity(y,pmax)

%% ADF Test
[hADFnotrend,pValnotrend,~,~,regADF] = adftest(y,'model','ARD','lags',0:pmax);
[hADFtrend,pValtrend,~,~,regADFtrend] = adftest(y,'model','TS','lags',0:pmax);

% Find optimal lag length for model without and with trend component
BICnotrend = [];
BICtrend = [];
for j=1:pmax+1
    BICnotrend = [BICnotrend,regADF(1,j).BIC];
end
for j=1:pmax+1
    BICtrend = [BICtrend,regADFtrend(1,j).BIC];
end
[BICminnotrend,pstarnotrend] = min(BICnotrend);
[BICmintrend,pstartrend] = min(BICtrend);

% Rejection decisions and p-values of ADF Test by lag length and deterministic terms
hADF = [hADFnotrend',hADFtrend'];
pValADF = [pValnotrend',pValtrend'];

% Rejection decisions and p-values of optimal models
hADF_notrend = hADF(pstarnotrend,1);
hADF_trend = hADF(pstartrend,2);
pValADF_notrend = pValADF(pstarnotrend,1);
pValADF_trend = pValADF(pstartrend,2);

%% KPSS Test
[hKPSSnotrend] = kpsstest(y,'Lags',0:pmax,'Trend',false);
[hKPSStrend] = kpsstest(y,'Lags',0:pmax,'Trend',true);
% Rejection decisions and p-values of KPSS Test by lag length and deterministic terms
hKPSS = [hKPSSnotrend',hKPSStrend'];
pValKPSS = [pValnotrend',pValtrend'];

% Rejection decisions and p-values of optimal models
hKPSS_notrend = hKPSSnotrend(pstarnotrend);
hKPSS_trend = hKPSStrend(pstartrend);
pValKPSS_notrend = pValKPSS(pstarnotrend,1);
pValKPSS_trend = pValKPSS(pstartrend,2);

%% Model comparison
if BICmintrend<BICminnotrend
    trend = 1;
else
    trend = 0;
end
trendpVal = regADFtrend(1,pstartrend).tStats.pVal(2);