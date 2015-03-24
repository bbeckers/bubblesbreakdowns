function [delta,BSADF,BSADFind,BSADFdates,cv_BSADF,tau0] = psy13(x,r0,adflag,mindur,bridge,startdate,cv)

%% Header
% This function performs the backward rolling windows sup ADF test of
% Philipps, Shi and Yu (2013) for detecting asset price bubbles and dating
% their emergence and collapse.
%_________________________________________________________________________%
% Inputs: x = Tx1 timeseries: asset prices or price-to-rent ratio in levels
%        r0 = Scalar [0,1]: ratio of the initial training sample to T
%    adflag = Scalar: lag length in the ADF regression equation                 
%    mindur = Scalar: number of periods a bubble must be present before the 
%             indicator switches to 1
%    bridge = Scalar: number of periods over which two otherwise seperate
%             bubbles are joined
% startdate = Matlab date: first observation of x
%        cv = Set of critical values adjusted for T and tau0
%_________________________________________________________________________%
% Outputs: BSADF = dim*1 timeseries: BSADF test statistic
%       BSADFind = dim*1 binary timeseries: 1 if bubble is present in
%                  period t, 0 otherwise
%     BSADFdates = Cell structure: Bubble emergence and collapse dates
%       cv_BSADF = dim*1 timeseries: BSADF critical values
%           tau0 = Scalar: Number of periods in initial training sample
%_________________________________________________________________________%

%% Function
if ~isnumeric(x)
    x=x(:);
end
x = log(x);

if cv==1
    % Simulate critical values
    simulation = 1;
    % Quantile(s) for critical values
    alpha = 0.95;
    % Number of simulations to obtain critical values
    sims = 2499;
else
    simulation = 0;
end

T = length(x);
% Initial training sample (first window length)
tau0 = floor(r0*T);
% End points
r2 = (tau0:1:T)';

% Size of monitoring period
dim = T-tau0+1;

% ADF-statistics
delta = zeros(dim,1);
BSADF = zeros(dim,1);
for t=1:dim
    % Current window lengths (minimal to maximal length)
    rwwindow = (tau0:1:r2(t))';
    % Current set of observations
    y = x(1:tau0-1+t);
    % Start points
    r1 = r2(t)-rwwindow+1;
    % ADF tests for all current windows
    b = zeros(size(rwwindow,1),1);
    adf_rw = zeros(size(rwwindow,1),1);
    for i=1:size(rwwindow,1)
        [b(i),~,~,adf_rw(i)] = ols_adf(y(r1(i):r2(t)),adflag);
    end
    % Collect surpremum of ADF tests for current windows
    [BSADF(t),ind] = max(adf_rw);
    delta(t) = b(ind);
end; clear i t

%% MC simulation of critical values
% CV BSADF
% If cv are not provided: Simulate
if simulation==1
    cv_BSADF = zeros(dim,1);
    for t=1:dim;
        t0 = r2(t);
        % Seed
        SI = 1;
        randn('seed',SI);
        % Generate random numbers for random walk with drift
        e = randn(t0,sims);
        a = t0^(-1);
        y = cumsum(e+a);
        % Test
        badfs = zeros(t0-tau0+1,sims);
        sadfs = ones(sims,1);
        for m=1:sims
            for i=tau0:1:t0
                [~,~,~,badfs(i-tau0+1,m)] = ols_adf(y(1:i,m),adflag);
            end
            sadfs(m,1) = max(badfs(:,m));
        end
        cv_BSADF(t) = quantile(sadfs,alpha);
    end;
    q = 2;
    cv_BSADF_MA = zeros(dim,1);
    for j=1:dim
        cv_BSADF_MA(j) = mean(cv_BSADF(max(j-q,1):min(j+q,dim)));
    end
else
    cv_BSADF = cv(1:dim);
end

%% Bubble indicator
BSADFind = (BSADF>cv_BSADF)*1;
% Bridging bubble periods
for j=1:dim
    if sum(BSADFind(j:min(j+bridge,dim)))>0 && BSADFind(max(1,j-1))>0
        BSADFind(j) = 1;
    end
end; clear j
% % Requiring bubbles to be at least [mindur] periods long
% BSADFind2 = zeros(dim,1);
% for j=max(mindur,1):1:dim
%     if mean(BSADFind(j-max(mindur,1)+1:j))==1
%         BSADFind2(j) = 1;
%     end
% end
% BSADFind = BSADFind2;

%% Switching dates
switchdates = zeros(length(BSADFind),1);
% Find observations where indicator variable switches
for i=2:length(BSADFind)
    if BSADFind(i)==BSADFind(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
% Convert index to date variables
index = find(switchdates==1);
BSADFdates = cell(length(index),1);
for i=1:sum(switchdates)
    BSADFdates{i} = dat2str(startdate+tau0+index(i)-2);
end; clear i index