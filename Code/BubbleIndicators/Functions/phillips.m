function [delta,ADF,cv,ADFind,ADFdates,tau0,deltaII] = phillips(x,r0,Jmax,bridge,H,grid,startdate,alpha,rolling,iie)
%% Header

%% Function
x = log(x);

if nargin<10
    iie = 0;
end
if nargin<9
    rolling = 0;
end
if nargin<8
    alpha = 0.05;
end

T = length(x);
tau0 = floor(r0*T);

delta = zeros(T-tau0+1,1);
if iie==1;
    deltaII = zeros(T-tau0+1,1);
end
ADF =  zeros(T-tau0+1,1);
cv = zeros(T-tau0+1,1);
beta = zeros(Jmax,T-tau0+1);

for i=1:T-tau0+1
    % Check for recursive or rolling estimation
    if rolling==1
        y = x(i:tau0-1+i);
    else y = x(1:tau0-1+i);
    end
    % Find optimal lag length by AIC or sequential testing
    AIC = zeros(Jmax,1);
    seqtest = zeros(Jmax,1);
    for j=1:Jmax
        [~,seqtest(j),AIC(j)] = ols(y,j,alpha);
    end; clear j
    seqtest_ind = find(seqtest,1,'first')-1;
    [~,AIC_ind] = min(AIC);
    Jstar = min(seqtest_ind,AIC_ind);
    % AR regression and ADF test
    [beta(1:Jstar+2,i),~,~,ADF(i)] = ols(y,Jstar,alpha);
    delta(i) = beta(2,i);
    % Indirect Inference estimator for delta
    if iie==1
        Omega = (round(grid^(-1)*delta(i))/grid^(-1):grid:delta(i)/0.97)';
        deltaII(i) = indirInf(y,beta(1:Jstar+2,i),Omega,H);
    end
    % Critical values of ADF test
    if rolling==1
        cv(i) = log(log(r0*T))/100;
    else cv(i) = log(log((r0+(i-1)/T)*T))/100;
    end
end; clear i Jstar check_ind AIC_ind AIC check

%% Bubble indicator
ADFind = (ADF>cv)*1;
% Bridging bubble periods that are less than [bridge] periods apart
for j=1:T-tau0+1
    if sum(ADFind(j:min(j+bridge-1,T-tau0+1)))>0 && ADFind(max(1,j-1))>0
        ADFind(j) = 1;
    end
end; clear j

%% Switching dates
switchdates = zeros(length(ADFind),1);
for i=2:length(ADFind)
    if ADFind(i)==ADFind(i-1)
        switchdates(i)=0;
    else switchdates(i)=1;
    end
end; clear i
index = find(switchdates==1);
ADFdates = cell(length(index),1);
for i=1:sum(switchdates)
    ADFdates{i} = dat2str(startdate+tau0+index(i)-2);
end; clear i index