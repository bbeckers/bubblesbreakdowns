function deltaII = indirInf(y,bLS,Omega,H)
%% Header
% This function computes the indirect inference estimate for the
% autoregressive parameter delta in the ADF-testing equation.
%
% Input:
% y: Tx1 time series
% bLS: Vector of parameters estimated by OLS from observations, deltaLS is
% second entry
% Omega: Nx1 vector of candidate autoregressive paramters from which
% deltaII is to be chosen
% H: Number of simulations for each delta in Omega
%
% Output:
% deltaII: scalar autoregressive parameter obtained from indirect inference

%% Function
T = length(y);
N = length(Omega);
if nargin<4
    H = 10000;
end

dist = zeros(N,1);
rng('default')
% Obtain H sequences of error terms each of length T
u = randn(T,H);
for n=1:N
    % Simulate H series of y for each possible value delta in Omega
    ysim = zeros(T,H);
    ysim(1:length(bLS)-1,:) = y(1:length(bLS)-1)*ones(1,H);
    deltasim = zeros(2,H);
    for h=1:H
        for j=length(bLS)+1:T
            dy = zeros(1,length(bLS)-2);
            for k=1:length(bLS)-2
                dy(k) = y(j-1)-y(j-1-k);
            end
            ysim(j,h) = bLS(1)+Omega(n)*ysim(j-1,h)+dy*bLS(3:end)+u(j,h);
        end
        % Estimate parameters from simulated time series
        X = [ones(T-1,1), ysim(1:T-1,h)];
        deltasim(:,h) = (X'*X)^(-1)*X'*ysim(2:T,h);
    end
    [~,msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:nearlySingularMatrix')
        dist(n:N) = 999;
        clear msgid; lastwarn('')
        break
    else deltasimbar = mean(deltasim,2);
        % Compute distance between delta estimated from true observations
        % and from simulated data
        dist(n) = sqrt((bLS(2)-deltasimbar(2))^2);
    end
end

% Match distance miniziming estimated delta and "true" data generating
% delta from set Omega
[~,ind] = min(dist);
deltaII = Omega(ind);