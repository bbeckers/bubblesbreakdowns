function [FLUC,FLUCind,FLUCdates,boundary,n] = fluc(x,r0,d,roll,mindur,bridge,startdate)
% -----------------------------------------------------------------------
% PURPOSE: performs fluctuation monitoring as in Homm and Breitung (2012)
%------------------------------------------------------------------------
% USAGE: rej = FLUC(y,n,k,d)
% where: x  = a time series vector of dimensions (T,1)
%        r0 = value between [0,1] determining size of initial training
%             sample in relation to T
%        d  = optional scalar indicating whether y shall be detrended
%             (d=1 in the case with detrending)
%        roll = optional scalar indicating whether recursive or rolling
%               estimation shall be employed (1 for rolling)
%        mindur = Minimum duration of ADF>cv for indicator jumping to 1
%        bridge = Scalar number of months that are bridged between two
%                 bubble periods
%        startdate = Date value of first observation for date stamping
%-----------------------------------------------------------------------
% RETURNS:
%        FLUC = (T-n)*1 binary timeseries: FLUC test statistic
%     FLUCind = (T-n)*1 binary timeseries: 1 if bubble is present in
%               period t, 0 otherwise
%   FLUCdates = Cell structure: Bubble emergence and collapse dates
%    boundary = (T-n)*1 timeseries: critical boundary
%           n = Scalar: Number of periods in initial training sample
%------------------------------------------------------------------------

x = log(x);

if nargin<4
    roll = 0;
end
if nargin<3
    d = 1;
end

T = length(x);

% Length of monitoring sample relative to training sample
k = 1/r0;
n = ceil(r0*T);

%% Critical values
% From Homm and Breitung (2012, Table 7)
if d==0 % Without detrending
    crit=[3.88 4.56 4.86 5.06 5.19 5.38 5.52
          4.19 4.80 5.11 5.34 5.50 5.72 5.81
          4.50 5.14 5.55 5.69 5.89 6.05 6.26];
else % With detrending
    crit=[5.50 6.52 7.24 7.96 8.04 8.45 9.13
          7.29 8.68 9.30 9.62 10.11 10.46 10.87
          8.12 9.82 10.45 10.82 11.20 11.55 11.80];
end
% Select critical value depending on length of training and monitoring
% sample
if k>=9
    j=7;
elseif k>=7 && k<9
    j=6;
else
    j=round(k-1);
end

if n <= 35
    i=1;
elseif n>35 && n<=75
    i=2;
else
    i=3;
end
const = crit(i,j);
% if r0==0.2
%     crit = 1.717008366435380;
% else
%     crit = 1.090957660493449;
% end
% const = crit;
boundary = sqrt(const+log((n+1:T)'/n));

%% ADF test
FLUC = zeros(T-n,1);
for t=1:T-n
    if roll==1
        y = x(t:n+t);
    else y = x(1:n+t);
    end
%     [~,~,~,FLUC(t)] = ols_adf(y,1);
    dimy = length(y);
    % Detrending by OLS
    if d==1
        X = [ones(dimy,1),(1:dimy)']; % Constant and time trend ((1:dimy)')
        a = X\y(1:dimy);
        z = y(1:dimy)-X*a;
    else
        z = y(1:dimy);
    end
    % Regressors and regressand
    dz = z(2:dimy)-z(1:dimy-1);
    z = z(1:dimy-1);
    % OLS regression
    b = z\dz;
    u = dz-z*b;
    sig2 = u'*u/(dimy-2);
    % ADF test statistic
    FLUC(t) = b*norm(z)/sqrt(sig2);
end
% Standardized ADF statistic (mean and standard deviation from Nabeya, 1999)
FLUC = (FLUC+2.1814)/0.7499;

%% Bubble indicator series
FLUCind = (FLUC>boundary)*1;
% Bridging bubble periods
for j=1:T-n
    if sum(FLUCind(j:min(j+bridge,T-n)))>0 && FLUCind(max(1,j-1))>0
        FLUCind(j) = 1;
    end
end; clear j
% % Account for persistence of bubble
% for j=1:T-n
%     FLUCind(T-j+1) = prod(FLUCind(max(1,T-n-j+2-par.mindur):T-n-j+1));
% end; clear j

%% Switching dates
switchdates = zeros(length(FLUCind),1);
for i=2:length(FLUCind)
    if FLUCind(i)==FLUCind(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
% Convert index to date variables
index = find(switchdates==1);
FLUCdates = cell(length(index),1);
for i=1:sum(switchdates)
    trash = dat2str(startdate+n+index(i)-1);
    FLUCdates{i} = trash{:};
end; clear i index