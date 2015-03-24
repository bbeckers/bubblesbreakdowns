function [delta,SADF,cv,SADFind,SADFdates,tau0,deltaII,Jstar] = pwy11(x,r0,Jmax,alpha,mindur,bridge,rolling,startdate,iie,H,grid)
%% Header
% This function performs the forward-recursive sup ADF test of
% Philipps, Shi and Yu (2011) for detecting asset price bubbles and dating
% their emergence and collapse.
%_________________________________________________________________________%
% Inputs: x = Tx1 timeseries: asset prices or price-to-rent ratio in levels
%        r0 = Scalar [0,1]: ratio of the initial training sample to T
%      Jmax = Scalar: maximum lag length in the ADF regression equation
%     alpha = Scalar: Significance level for sequential testing of 
%             estimated AR-model                 
%    mindur = Scalar: number of periods a bubble must be present before the 
%             indicator switches to 1
%    bridge = Scalar: number of periods over which two otherwise seperate
%             bubbles are joined
%   rolling = Binary variable: 1 if rolling windows shall be employed
% startdate = Matlab date: first observation of x
%       iie = Binary variable: 1 if indirect inference estimator (IIE)
%             shall be obtained
%         H = Scalar: Number of simulations for IIE
%      grid = Scalar: Grid point density for IIE (distance between points)
%_________________________________________________________________________%
% Outputs: delta = dim*1 timeseries: Estimated AR-parameter
%           SADF = dim*1 binary timeseries: Sup ADF test statistic
%             cv = dim*1 timeseries: SADF critical values
%        SADFind = dim*1 binary timeseries: 1 if bubble is present in
%                  period t, 0 otherwise
%     ´SADFdates = Cell structure: Bubble emergence and collapse dates
%           tau0 = Scalar: Number of periods in initial training sample
%        deltaII = dim*1 timeseries: Estimated AR-parameter by indirect
%                  inference
%          Jstar = dim*1 timeseries: Optimal lag length in period t
%_________________________________________________________________________%

%% Function
if ~isnumeric(x)
    x=x(:);
end
x = log(x);

if nargin<9
    iie = 0;
end

T = length(x);
tau0 = floor(r0*T);

% Size of monitoring period
dim = T-tau0+1;

delta = zeros(dim,1);
deltaII = zeros(dim,1);
SADF =  zeros(dim,1);
cv = zeros(dim,1);
beta = zeros(Jmax,dim);
Jstar = zeros(dim,1);

for i=1:dim
    % Check for recursive or rolling estimation
    if rolling==1
        y = x(i:tau0-1+i);
    else
        y = x(1:tau0-1+i);
    end
    % Find optimal lag length by AIC or sequential testing
    if Jmax>0
        AIC = zeros(Jmax,1);
        seqtest = zeros(Jmax,1);
        for j=1:Jmax
            [~,seqtest(j),AIC(j)] = ols_adf(y,j,alpha);
        end; clear j
        if any(seqtest)>0;
            seqtest_ind = find(seqtest,1,'last');
        else
            seqtest_ind = 0;
        end
        [~,AIC_ind] = min(AIC);
        Jstar(i) = AIC_ind;
    else
        Jstar(i) = 0;
    end
    % AR regression and ADF test with optimal lag length
    [beta(1:Jstar(i)+2,i),~,~,SADF(i)] = ols_adf(y,Jstar(i),alpha);
    delta(i) = beta(2,i);
    % Indirect Inference estimator for delta
    if iie==1
        Omega = (round(grid^(-1)*delta(i))/grid^(-1):grid:delta(i)/0.97)';
        deltaII(i) = indirInf(y,beta(1:Jstar(i)+2,i),Omega,H);
    else
        deltaII(i) = delta(i);
    end
    % Critical values of ADF test
    if rolling==1
        cv(i) = log(log(r0*T))/100;
    else cv(i) = log(log((r0+(i-1)/T)*T))/100;
    end
end; clear i

%% Bubble indicator
SADFind = (SADF>cv)*1;
% Bridging bubble periods that are less than [bridge] periods apart
for j=1:dim
    if sum(SADFind(j:min(j+bridge,dim)))>0 && SADFind(max(1,j-1))>0
        SADFind(j) = 1;
    end
end; clear j
% % Account for persistence of bubble
% for j=1:T
%     SADFind(T-j+1) = prod(SADFind(max(1,T-j+2-par.mindur):T-j+1));
% end; clear j

%% Switching dates
switchdates = zeros(length(SADFind),1);
for i=2:length(SADFind)
    if SADFind(i)==SADFind(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
% Convert index to date variables
index = find(switchdates==1);
SADFdates = cell(length(index),1);
for i=1:sum(switchdates)
    SADFdates{i} = dat2str(startdate+tau0+index(i)-2);
end; clear i index