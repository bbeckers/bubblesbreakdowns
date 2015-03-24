%% Header
clear
close all
clc
addpath(genpath('Functions'));
addpath(genpath('IRIS_Tbx'));
% irisstartup

%% Specifications
assetname = 'Stock Prices'; % Choose between 'Stock Prices' and 'House Prices'
if strcmp(assetname,'Stock Prices')
    asset = 's';
else asset = 'h';
end

%% Load Data
D = dbload('..\Daten\Non-Revised\NonrevData.csv','nameRow=',1);

%% Date and sampling
par = struct();
par.freq = D.CPI.freq;
rangeunb = dbrange(D,fieldnames(D),'startdate','unbalanced','enddate','unbalanced');
date = struct();
if strcmp(asset,'s')
    date.start = rangeunb(1);
    date.end = rangeunb(end);
else
    date.start = get(D.DISPINC,'start');
    date.end = get(D.HPI,'end');
end
% Collect dates as strings
clear rangeunb
date.start = mm(1959,1);
date.end = mm(2014,08);
a = dat2str(date.start);
date.starty = a{1,1}(1:4); clear a
date.starty  = str2double(date.starty);
T = date.end-date.start+1;
timeline = (date.starty+1/par.freq:1/par.freq:date.starty+T/par.freq)';

%% Data transformation
% Set base year to 2000M1
D.CPI = 100*D.CPI/mean(D.CPI(mm(2013,1):mm(2013,12)));
% D.SP500 = 100*D.SP500/mean(D.SP500(mm(2013,1):mm(2013,12)));
% D.DIV = 100*D.DIV/mean(D.DIV(mm(2013,1):mm(2013,12)));
D.EAR = 100*D.EAR/mean(D.EAR(mm(2013,1):mm(2013,12)));
% D.HPI = 100*D.HPI/mean(D.HPI(mm(2013,1):mm(2013,12)));
% D.RENT = 100*D.RENT/mean(D.RENT(mm(2013,1):mm(2013,12)));
% D.RDISPINCCAP = 100*D.RDISPINCCAP/mean(D.RDISPINCCAP(mm(2013,1):mm(2013,12)));

% Compute real stock prices, dividends and earnings
if strcmp(asset,'s')
    p = 100*D.SP500./D.CPI;
    d = 100*D.DIV./D.CPI;
    e = 100*D.EAR./D.CPI;
else
    p = 100*D.HPI./D.CPI;
    d = 100*D.RENT./D.CPI;
    e = D.RDISPINCCAP;
end

p = tseries(date.start:date.end,p(date.start:date.end));

if strcmp(asset,'s')
    f = tseries(date.start:date.end,d(date.start:date.end));
else
    f = tseries(date.start:date.end,e(date.start:date.end));
end
p2f = p./f;

%% Parametrization
% Phillips 2011
% Maximum lag order in autoregressive equation for Phillips-Test
par.Jmax = 12;
% Significance level for sequential testing of estimated AR-model (ADF)
par.alpha.ols = 0.05;
% Compute indirect inference estimator for autoregressive parameter? (0 for
% no, 1 for yes)
par.iie = 0;
% Density of grid over which ii estimator shall be evaluated
par.grid = 0.001;
% Indirect inference simulations (#)
par.H = 1000;

% Phillips 2011 - rolling
% Rolling or recursive estimation of AR (1 for rolling)
par.rollingl = 1;
% Initial training period for estimation
par.r0l = 0.2;

% Phillips 2011 - recursive
% Rolling or recursive estimation of AR (1 for rolling)
par.rollings = 0;
% Initial training period for estimation
par.r0s = 0.01 + 1.8/sqrt(T);

% Phillips 2013
% Lag order in autoregressive equation for Phillips-Test (FIXED)
par.adflag = 1;
% Initial training period for estimation
par.r0 = 0.01+1.8/sqrt(T);

% FLUC Test
% Detrending
par.detrend = 1;

% HP-filter
% Smoothing parameter
par.lambda = 100000*(par.freq/4)^4;
% Bubble threshold
if strcmp(asset,'s')
    par.hpthreshold = 0.1;
else
    par.hpthreshold = 0.075;
end
% Rolling window length
par.w = 8*par.freq;

% Minimum duration of BSADF>cv for bubble indication (in months)
par.mindur = 0;%floor(log(T));
% Bridging period
par.bridge = floor(log(T));

%% Asset price bubbles - Testing for explosiveness

%% Phillips 2011 - recursive
[~,PWY11sp,cvPWY11s,IndPWY11sp,DatesPWY11sp,tau0PWY11s] = pwy11(p,par.r0s,par.Jmax,par.alpha.ols,par.mindur,0,par.rollings,date.start,par.iie,par.H,par.grid);
[~,PWY11sf,~,IndPWY11sf,DatesPWY11sf] = pwy11(f,par.r0s,par.Jmax,par.alpha.ols,par.mindur,0,par.rollings,date.start,par.iie,par.H,par.grid);
% Keep price bubbles only
IndPWY11sindiv = IndPWY11sp-IndPWY11sf;
IndPWY11sindiv(IndPWY11sindiv==-1) = 0;
% Bridging bubble periods
for j=1:length(IndPWY11sindiv)
    if sum(IndPWY11sindiv(j:min(j+par.bridge,length(IndPWY11sindiv))))>0 && IndPWY11sindiv(max(1,j-1))>0
        IndPWY11sindiv(j) = 1;
    end
end; clear j
% Find dates
switchdates = IndPWY11sindiv(2:end)-IndPWY11sindiv(1:end-1);
switchdates = find(switchdates~=0);
DatesPWY11s = dat2str(date.start+tau0PWY11s-1+switchdates);
[~,PWY11sr,~,IndPWY11sratio,DatesPWY11sr] = pwy11(p2f,par.r0s,par.Jmax,par.alpha.ols,par.mindur,par.bridge,par.rollings,date.start,par.iie,par.H,par.grid);
% Combine bubbles in individual series and in ratio
IndPWY11scomb = ((IndPWY11sindiv+IndPWY11sratio)>0)*1;
switchdates = IndPWY11scomb(2:end)-IndPWY11scomb(1:end-1);
switchdates = find(switchdates~=0);
DatesPWY11scomb = dat2str(date.start+tau0PWY11s-1+switchdates);

IndPWY11sp = [NaN(tau0PWY11s-1,1);IndPWY11sp];
IndPWY11sf = [NaN(tau0PWY11s-1,1);IndPWY11sf];
IndPWY11sindiv = [NaN(tau0PWY11s-1,1);IndPWY11sindiv];
IndPWY11sratio = [NaN(tau0PWY11s-1,1);IndPWY11sratio];
IndPWY11scomb = [NaN(tau0PWY11s-1,1);IndPWY11scomb];

%% Phillips 2011 - rolling
[~,PWY11lp,cvPWY11l,IndPWY11lp,DatesPWY11lp,tau0PWY11l] = pwy11(p,par.r0l,par.Jmax,par.alpha.ols,par.mindur,0,par.rollingl,date.start,par.iie,par.H,par.grid);
[~,PWY11lf,~,IndPWY11lf,DatesPWY11lf] = pwy11(f,par.r0l,par.Jmax,par.alpha.ols,par.mindur,0,par.rollingl,date.start,par.iie,par.H,par.grid);
% Keep price bubbles only
IndPWY11lindiv = IndPWY11lp-IndPWY11lf;
IndPWY11lindiv(IndPWY11lindiv==-1) = 0;
% Bridging bubble periods
for j=1:length(IndPWY11lindiv)
    if sum(IndPWY11lindiv(j:min(j+par.bridge,length(IndPWY11lindiv))))>0 && IndPWY11lindiv(max(1,j-1))>0
        IndPWY11lindiv(j) = 1;
    end
end; clear j
% Find dates
switchdates = IndPWY11lindiv(2:end)-IndPWY11lindiv(1:end-1);
switchdates = find(switchdates~=0);
DatesPWY11l = dat2str(date.start+tau0PWY11l-1+switchdates);
[~,PWY11lr,~,IndPWY11lratio,DatesPWY11lr] = pwy11(p2f,par.r0l,par.Jmax,par.alpha.ols,par.mindur,par.bridge,par.rollingl,date.start,par.iie,par.H,par.grid);
% Combine bubbles in individual series and in ratio
IndPWY11lcomb = ((IndPWY11lindiv+IndPWY11lratio)>0)*1;
switchdates = IndPWY11lcomb(2:end)-IndPWY11lcomb(1:end-1);
switchdates = find(switchdates~=0);
DatesPWY11lcomb = dat2str(date.start+tau0PWY11l-1+switchdates);

IndPWY11lp = [NaN(tau0PWY11l-1,1);IndPWY11lp];
IndPWY11lf = [NaN(tau0PWY11l-1,1);IndPWY11lf];
IndPWY11lindiv = [NaN(tau0PWY11l-1,1);IndPWY11lindiv];
IndPWY11lratio = [NaN(tau0PWY11l-1,1);IndPWY11lratio];
IndPWY11lcomb = [NaN(tau0PWY11l-1,1);IndPWY11lcomb];

%% Phillips 2013
% Load critical values
load('Functions\Shi\cv_bsadf_670.mat','quantile_badfs_MA');
cv_bsadf = quantile_badfs_MA; clear quantile_badfs_MA
[~,PSY13p,IndPSY13p,DatesPSY13p,~,tau0PSY13] = psy13(p,par.r0,par.adflag,par.mindur,0,date.start,cv_bsadf);
[~,PSY13f,IndPSY13f,DatesPSY13f] = psy13(f,par.r0,par.adflag,par.mindur,0,date.start,cv_bsadf);
% Keep price bubbles only
IndPSY13indiv = IndPSY13p-IndPSY13f;
IndPSY13indiv(IndPSY13indiv==-1) = 0;
% Bridging bubble periods
for j=1:length(IndPSY13indiv)
    if sum(IndPSY13indiv(j:min(j+par.bridge,length(IndPSY13indiv))))>0 && IndPSY13indiv(max(1,j-1))>0
        IndPSY13indiv(j) = 1;
    end
end; clear j
% Find dates
switchdates = IndPSY13indiv(2:end)-IndPSY13indiv(1:end-1);
switchdates = find(switchdates~=0);
DatesPSY13 = dat2str(date.start+tau0PSY13-1+switchdates);
[~,PSY13r,IndPSY13ratio,DatesPSY13r] = psy13(p2f,par.r0,par.adflag,par.mindur,par.bridge,date.start,cv_bsadf);
% Combine bubbles in individual series and in ratio
IndPSY13comb = ((IndPSY13indiv+IndPSY13ratio)>0)*1;
switchdates = IndPSY13comb(2:end)-IndPSY13comb(1:end-1);
switchdates = find(switchdates~=0);
DatesPSY13comb = dat2str(date.start+tau0PSY13-1+switchdates);

IndPSY13p = [NaN(tau0PSY13-1,1);IndPSY13p];
IndPSY13f = [NaN(tau0PSY13-1,1);IndPSY13f];
IndPSY13indiv = [NaN(tau0PSY13-1,1);IndPSY13indiv];
IndPSY13ratio = [NaN(tau0PSY13-1,1);IndPSY13ratio];
IndPSY13comb = [NaN(tau0PSY13-1,1);IndPSY13comb];

%% FLUC Test
par.r0s = 0.1;

% Recursive, short
[FLUCsp,IndFLUCsp,DatesFLUCsp,~,tau0FLUCs] = fluc(p(:),par.r0s,par.detrend,0,par.mindur,0,date.start);
[FLUCsf,IndFLUCsf,DatesFLUCsf] = fluc(f(:),par.r0s,par.detrend,0,par.mindur,0,date.start);
% Keep price bubbles only
IndFLUCsindiv = IndFLUCsp-IndFLUCsf;
IndFLUCsindiv(IndFLUCsindiv==-1) = 0;
% Bridging bubble periods
for j=1:length(IndFLUCsindiv)
    if sum(IndFLUCsindiv(j:min(j+par.bridge,length(IndFLUCsindiv))))>0 && IndFLUCsindiv(max(1,j-1))>0
        IndFLUCsindiv(j) = 1;
    end
end; clear j
% Find dates
switchdates = IndFLUCsindiv(2:end)-IndFLUCsindiv(1:end-1);
switchdates = find(switchdates~=0);
DatesFLUCs = dat2str(date.start+tau0FLUCs-2+switchdates);
[FLUCsr,IndFLUCsratio,DatesFLUCsr,] = fluc(p2f(:),par.r0s,par.detrend,0,par.mindur,par.bridge,date.start);
% Combine bubbles in individual series and in ratio
IndFLUCscomb = ((IndFLUCsindiv+IndFLUCsratio)>0)*1;
switchdates = IndFLUCscomb(2:end)-IndFLUCscomb(1:end-1);
switchdates = find(switchdates~=0);
DatesFLUCscomb = dat2str(date.start+tau0FLUCs-1+switchdates);

IndFLUCsp = [NaN(tau0FLUCs,1);IndFLUCsp];
IndFLUCsf = [NaN(tau0FLUCs,1);IndFLUCsf];
IndFLUCsindiv = [NaN(tau0FLUCs,1);IndFLUCsindiv];
IndFLUCsratio = [NaN(tau0FLUCs,1);IndFLUCsratio];
IndFLUCscomb = [NaN(tau0FLUCs,1);IndFLUCscomb];

% Recursive, long
[FLUClp,IndFLUClp,DatesFLUClp,~,tau0FLUCl] = fluc(p(:),par.r0s,par.detrend,0,par.mindur,0,date.start);
[FLUClf,IndFLUClf,DatesFLUClf] = fluc(f(:),par.r0s,par.detrend,0,par.mindur,0,date.start);
% Keep price bubbles only
IndFLUClindiv = IndFLUClp-IndFLUClf;
IndFLUClindiv(IndFLUClindiv==-1) = 0;
% Bridging bubble periods
for j=1:length(IndFLUClindiv)
    if sum(IndFLUClindiv(j:min(j+par.bridge,length(IndFLUClindiv))))>0 && IndFLUClindiv(max(1,j-1))>0
        IndFLUClindiv(j) = 1;
    end
end; clear j
% Find dates
switchdates = IndFLUClindiv(2:end)-IndFLUClindiv(1:end-1);
switchdates = find(switchdates~=0);
DatesFLUCl = dat2str(date.start+tau0FLUCl-2+switchdates);
[FLUClr,IndFLUClratio,DatesFLUClr,] = fluc(p2f(:),par.r0s,par.detrend,0,par.mindur,par.bridge,date.start);
% Combine bubbles in individual series and in ratio
IndFLUClcomb = ((IndFLUClindiv+IndFLUClratio)>0)*1;
switchdates = IndFLUClcomb(2:end)-IndFLUClcomb(1:end-1);
switchdates = find(switchdates~=0);
DatesFLUClcomb = dat2str(date.start+tau0FLUCl-1+switchdates);

IndFLUClp = [NaN(tau0FLUCs,1);IndFLUClp];
IndFLUClf = [NaN(tau0FLUCl,1);IndFLUClf];
IndFLUClindiv = [NaN(tau0FLUCl,1);IndFLUClindiv];
IndFLUClratio = [NaN(tau0FLUCl,1);IndFLUClratio];
IndFLUClcomb = [NaN(tau0FLUCl,1);IndFLUClcomb];

%% HP-Filter
par.mindur = 1;

% Recursive
ptrend = one_sided_hp_filter_serial(log(p(:)),par.lambda);
HPrec = log(p(:))-ptrend; % Size of Bubble
IndHPrec = (HPrec>par.hpthreshold)*1;
% Bridging bubble periods
for j=1:T
    if sum(IndHPrec(j:min(j+par.bridge,T)))>0 && IndHPrec(max(1,j-1))>0
        IndHPrec(j) = 1;
    end
end; clear j
% Account for persistence of bubble
for j=1:T
    IndHPrec(T-j+1) = prod(IndHPrec(max(1,T-j+2-par.mindur):T-j+1));
end; clear j
% Switching dates
switchdates = zeros(length(IndHPrec),1);
for i=2:length(IndHPrec)
    if IndHPrec(i)==IndHPrec(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
% Convert index to date variables
index = find(switchdates==1);
DatesHPrec = cell(length(index),1);
for i=1:sum(switchdates)
    trash = dat2str(date.start+index(i)-1);
    DatesHPrec{i} = trash{:};
end; clear i index trash

% Rolling
ptrend = zeros(T,1);
ptrend(1:par.w) = one_sided_hp_filter_serial(log(p(date.start:date.start+par.w-1)),par.lambda);
for j=par.w+1:1:T
    ptrendcur = one_sided_hp_filter_serial(log(p(date.start+j-par.w:date.start+j-1)),par.lambda);
    ptrend(j) = ptrendcur(end);
end
HProl = log(p(:))-ptrend; % Size of Bubble
IndHProl = (HProl>par.hpthreshold)*1;
% Bridging bubble periods
for j=1:T
    if sum(IndHProl(j:min(j+par.bridge,T)))>0 && IndHProl(max(1,j-1))>0
        IndHProl(j) = 1;
    end
end; clear j
% Account for persistence of bubble
for j=1:T
    IndHProl(T-j+1) = prod(IndHProl(max(1,T-j+2-par.mindur):T-j+1));
end; clear j
% Switching dates
switchdates = zeros(length(IndHProl),1);
for i=2:length(IndHProl)
    if IndHProl(i)==IndHProl(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
% Convert index to date variables
index = find(switchdates==1);
DatesHProl = cell(length(index),1);
for i=1:sum(switchdates)
    DatesHProl{i} = dat2str(date.start+index(i)-1);
end; clear i index

%% Combinations
% All indicators
IndComb1 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=1)*1;
IndComb2 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=2)*1;
IndComb3 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=3)*1;
IndComb4 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=4)*1;
IndComb5 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=5)*1;
IndComb6 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=6)*1;
IndComb7 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=7)*1;
IndComb8 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio+IndFLUCsindiv+IndFLUCsratio+IndFLUClindiv+IndFLUClratio+IndHPrec+IndHProl)>=8)*1;

% Phillips only
IndCombPhil1 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio)>=1)*1;
IndCombPhil2 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio)>=2)*1;
IndCombPhil3 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio)>=3)*1;
IndCombPhil4 = ((IndPWY11sindiv+IndPWY11sratio+IndPWY11lindiv+IndPWY11lratio+IndPSY13indiv+IndPSY13ratio)>=4)*1;

% Best detectors
IndCombBest1 = ((IndPWY11sindiv+IndPWY11lindiv+IndPSY13indiv+IndHPrec)>=1)*1;
IndCombBest2 = ((IndPWY11sindiv+IndPWY11lindiv+IndPSY13indiv+IndHPrec)>=2)*1;
IndCombBest3 = ((IndPWY11sindiv+IndPWY11lindiv+IndPSY13indiv+IndHPrec)>=3)*1;
IndCombBest4 = ((IndPWY11sindiv+IndPWY11lindiv+IndPSY13indiv+IndHPrec)>=4)*1;

switchdates1 = IndComb1(2:end)-IndComb1(1:end-1);
switchdates1 = find(switchdates1~=0);
DatesComb1 = dat2str(date.start+switchdates1);
switchdates2 = IndComb2(2:end)-IndComb2(1:end-1);
switchdates2 = find(switchdates2~=0);
DatesComb2 = dat2str(date.start+switchdates2);
switchdates3 = IndComb3(2:end)-IndComb3(1:end-1);
switchdates3 = find(switchdates3~=0);
DatesComb3 = dat2str(date.start+switchdates3);
switchdates4 = IndComb4(2:end)-IndComb4(1:end-1);
switchdates4 = find(switchdates4~=0);
DatesComb4 = dat2str(date.start+switchdates4);
switchdates5 = IndComb5(2:end)-IndComb5(1:end-1);
switchdates5 = find(switchdates5~=0);
DatesComb5 = dat2str(date.start+switchdates5);
switchdates6 = IndComb6(2:end)-IndComb6(1:end-1);
switchdates6 = find(switchdates6~=0);
DatesComb6 = dat2str(date.start+switchdates6);
switchdates7 = IndComb6(2:end)-IndComb6(1:end-1);
switchdates7 = find(switchdates7~=0);
DatesComb7 = dat2str(date.start+switchdates7);
switchdates8 = IndComb6(2:end)-IndComb6(1:end-1);
switchdates8 = find(switchdates8~=0);
DatesComb8 = dat2str(date.start+switchdates8);

switchdatesPhil1 = IndCombPhil1(2:end)-IndCombPhil1(1:end-1);
switchdatesPhil1 = find(switchdatesPhil1~=0);
DatesCombPhil1 = dat2str(date.start+switchdatesPhil1);
switchdatesPhil2 = IndCombPhil2(2:end)-IndCombPhil2(1:end-1);
switchdatesPhil2 = find(switchdatesPhil2~=0);
DatesCombPhil2 = dat2str(date.start+switchdatesPhil2);
switchdatesPhil3 = IndCombPhil3(2:end)-IndCombPhil3(1:end-1);
switchdatesPhil3 = find(switchdatesPhil3~=0);
DatesCombPhil3 = dat2str(date.start+switchdatesPhil3);
switchdatesPhil4 = IndCombPhil4(2:end)-IndCombPhil4(1:end-1);
switchdatesPhil4 = find(switchdatesPhil4~=0);
DatesCombPhil4 = dat2str(date.start+switchdatesPhil4);

switchdatesBest1 = IndCombBest1(2:end)-IndCombBest1(1:end-1);
switchdatesBest1 = find(switchdatesBest1~=0);
DatesCombBest1 = dat2str(date.start+switchdatesBest1);
switchdatesBest2 = IndCombBest2(2:end)-IndCombBest2(1:end-1);
switchdatesBest2 = find(switchdatesBest2~=0);
DatesCombBest2 = dat2str(date.start+switchdatesBest2);
switchdatesBest3 = IndCombBest3(2:end)-IndCombBest3(1:end-1);
switchdatesBest3 = find(switchdatesBest3~=0);
DatesCombBest3 = dat2str(date.start+switchdatesBest3);
switchdatesBest4 = IndCombBest4(2:end)-IndCombBest4(1:end-1);
switchdatesBest4 = find(switchdatesBest4~=0);
DatesCombBest4 = dat2str(date.start+switchdatesBest4);

%% Convert to timeseries
IndPWY11sp = tseries(date.start:date.end,IndPWY11sp);
IndPWY11sf = tseries(date.start:date.end,IndPWY11sf);
IndPWY11sindiv = tseries(date.start:date.end,IndPWY11sindiv);
IndPWY11sratio = tseries(date.start:date.end,IndPWY11sratio);
IndPWY11scomb = tseries(date.start:date.end,IndPWY11scomb);

IndPWY11lp = tseries(date.start:date.end,IndPWY11lp);
IndPWY11lf = tseries(date.start:date.end,IndPWY11lf);
IndPWY11lindiv = tseries(date.start:date.end,IndPWY11lindiv);
IndPWY11lratio = tseries(date.start:date.end,IndPWY11lratio);
IndPWY11lcomb = tseries(date.start:date.end,IndPWY11lcomb);

IndPSY13p = tseries(date.start:date.end,IndPSY13p);
IndPSY13f = tseries(date.start:date.end,IndPSY13f);
IndPSY13indiv = tseries(date.start:date.end,IndPSY13indiv);
IndPSY13ratio = tseries(date.start:date.end,IndPSY13ratio);
IndPSY13comb = tseries(date.start:date.end,IndPSY13comb);

IndFLUCsp = tseries(date.start:date.end,IndFLUCsp);
IndFLUCsf = tseries(date.start:date.end,IndFLUCsf);
IndFLUCsindiv = tseries(date.start:date.end,IndFLUCsindiv);
IndFLUCsratio = tseries(date.start:date.end,IndFLUCsratio);
IndFLUCscomb = tseries(date.start:date.end,IndFLUCscomb);

IndFLUClp = tseries(date.start:date.end,IndFLUClp);
IndFLUClf = tseries(date.start:date.end,IndFLUClf);
IndFLUClindiv = tseries(date.start:date.end,IndFLUClindiv);
IndFLUClratio = tseries(date.start:date.end,IndFLUClratio);
IndFLUClcomb = tseries(date.start:date.end,IndFLUClcomb);

IndHPrec = tseries(date.start:date.end,IndHPrec);
IndHProl = tseries(date.start:date.end,IndHProl);

IndComb1 = tseries(date.start:date.end,IndComb1);
IndComb2 = tseries(date.start:date.end,IndComb2);
IndComb3 = tseries(date.start:date.end,IndComb3);
IndComb4 = tseries(date.start:date.end,IndComb4);
IndComb5 = tseries(date.start:date.end,IndComb5);
IndComb6 = tseries(date.start:date.end,IndComb6);
IndComb7 = tseries(date.start:date.end,IndComb7);
IndComb8 = tseries(date.start:date.end,IndComb8);

IndCombPhil1 = tseries(date.start:date.end,IndCombPhil1);
IndCombPhil2 = tseries(date.start:date.end,IndCombPhil2);
IndCombPhil3 = tseries(date.start:date.end,IndCombPhil3);
IndCombPhil4 = tseries(date.start:date.end,IndCombPhil4);

IndCombBest1 = tseries(date.start:date.end,IndCombBest1);
IndCombBest2 = tseries(date.start:date.end,IndCombBest2);
IndCombBest3 = tseries(date.start:date.end,IndCombBest3);
IndCombBest4 = tseries(date.start:date.end,IndCombBest4);

%% Results
clear -regexp switch
% Figures

figure(1)
bar(timeline,IndCombBest2(:)*max(10000*p2f),'FaceColor',[0.5,0.5,0.5])
hold on
plot(timeline,10000*p2f(:),'LineWidth',1.5)
axis tight
title('House price bubble periods','FontSize',14)
ylabel('Price-to-Income Ratio','FontSize',12)

% Excel-Output
R = struct();
list = who('Ind*');
for i=1:length(list)
    R.(list{i}) = eval(list{i});
end

filename = strcat('..\Results\Indicators\',assetname,'\Phillips.csv');
dbsave(R,filename);
clear filename

save(strcat('..\Results\Indicators\',assetname,'\PSY13indiv.mat'),'PSY13p','PSY13f','IndPSY13indiv','DatesPSY13','tau0PSY13','p','f','cv_bsadf')
save(strcat('..\Results\Indicators\',assetname,'\PSY13ratio.mat'),'PSY13r','IndPSY13ratio','DatesPSY13r','tau0PSY13','p','f','cv_bsadf')
save(strcat('..\Results\Indicators\',assetname,'\PSY13comb.mat'),'IndPSY13comb','DatesPSY13comb','tau0PSY13','p','f')
save(strcat('..\Results\Indicators\',assetname,'\PWY11lindiv.mat'),'PWY11lp','PWY11lf','IndPWY11lindiv','DatesPWY11l','tau0PWY11l','p','f','cvPWY11l')
save(strcat('..\Results\Indicators\',assetname,'\PWY11lratio.mat'),'PWY11lr','IndPWY11lratio','DatesPWY11lr','tau0PWY11l','p','f','cvPWY11l')
save(strcat('..\Results\Indicators\',assetname,'\PWY11lcomb.mat'),'IndPWY11lcomb','DatesPWY11lcomb','tau0PWY11l','p','f')
save(strcat('..\Results\Indicators\',assetname,'\PWY11sindiv.mat'),'PWY11sp','PWY11sf','IndPWY11sindiv','DatesPWY11s','tau0PWY11s','p','f','cvPWY11s')
save(strcat('..\Results\Indicators\',assetname,'\PWY11sratio.mat'),'PWY11sr','IndPWY11sratio','DatesPWY11sr','tau0PWY11s','p','f','cvPWY11s')
save(strcat('..\Results\Indicators\',assetname,'\PWY11scomb.mat'),'IndPWY11scomb','DatesPWY11scomb','tau0PWY11s','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUCsindiv.mat'),'FLUCsp','FLUCsf','IndFLUCsindiv','DatesFLUCs','tau0FLUCs','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUCsratio.mat'),'FLUCsr','IndFLUCsratio','DatesFLUCsr','tau0FLUCs','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUCscomb.mat'),'IndFLUCscomb','DatesFLUCscomb','tau0FLUCs','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUClindiv.mat'),'FLUClp','FLUClf','IndFLUClindiv','DatesFLUCl','tau0FLUCl','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUClratio.mat'),'FLUClr','IndFLUClratio','DatesFLUClr','tau0FLUCl','p','f')
save(strcat('..\Results\Indicators\',assetname,'\FLUClcomb.mat'),'IndFLUClcomb','DatesFLUClcomb','tau0FLUCl','p','f')
save(strcat('..\Results\Indicators\',assetname,'\HPrec.mat'),'HPrec','IndHPrec','DatesHPrec')
save(strcat('..\Results\Indicators\',assetname,'\HProl.mat'),'HProl','IndHProl','DatesHProl')
save(strcat('..\Results\Indicators\',assetname,'\Comb1.mat'),'IndComb1','DatesComb1','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb2.mat'),'IndComb2','DatesComb2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb3.mat'),'IndComb3','DatesComb3','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb4.mat'),'IndComb4','DatesComb4','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb5.mat'),'IndComb5','DatesComb5','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb6.mat'),'IndComb6','DatesComb6','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb7.mat'),'IndComb7','DatesComb7','p','f')
save(strcat('..\Results\Indicators\',assetname,'\Comb8.mat'),'IndComb8','DatesComb8','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombPhil1.mat'),'IndCombPhil1','DatesCombPhil1','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombPhil2.mat'),'IndCombPhil2','DatesCombPhil2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombPhil3.mat'),'IndCombPhil3','DatesCombPhil2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombPhil4.mat'),'IndCombPhil4','DatesCombPhil2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombBest1.mat'),'IndCombBest1','DatesCombBest1','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombBest2.mat'),'IndCombBest2','DatesCombBest2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombBest3.mat'),'IndCombBest3','DatesCombBest2','p','f')
save(strcat('..\Results\Indicators\',assetname,'\CombBest4.mat'),'IndCombBest4','DatesCombBest2','p','f')