function [sample,lp_Sample,len] ...
    = regen2(This,NDraw,varargin)
% arwm  Regeneration time MCMC Metropolis posterior simulator.
%
% Syntax
% =======
%
%     [Theta,LogPost,AR,Scale,FinalCov] = regen(Pos,NDraw,...)
%
% Input arguments
% ================
%
% * `Pos` [ poster ] - Initialised posterior simulator object.
%
% * `NDraw` [ numeric ] - Length of the chain not including burn-in.
%
% Output arguments
% =================
%
% * `Theta` [ numeric ] - MCMC chain with individual parameters in rows.
%
% * `LogPost` [ numeric ] - Vector of log posterior density (up to a
% constant) in each draw.
%
% * `AR` [ numeric ] - Vector of cumulative acceptance ratios in each draw.
%
% * `Scale` [ numeric ] - Vector of proposal scale factors in each draw.
%
% * `FinalCov` [ numeric ] - Final proposal covariance matrix; the final
% covariance matrix of the random walk step is Scale(end)^2*FinalCov.
%
% Options
% ========
%
% References
% ========
% 1. Brockwell, A.E., and Kadane, J.B., 2004. "Identification of ]
%    Regeneration Times in MCMC Simulation, with Application to Adaptive
%    Schemes," mimeo, Carnegie Mellon University.
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team & Bojan Bejanov & Troy Matheson.

% Validate required inputs.
pp = inputParser();
pp.addRequired('Pos',@(x) isa(x,'poster'));
pp.addRequired('NDraw',@isnumericscalar);
pp.parse(This,NDraw);

% Parse options.
opt = passvalopt('poster.regen',varargin{:});

%--------------------------------------------------------------------------

s = mylogpoststruct(This);

if opt.initialChainSize < 1
    % initial chain size is a percentage
    opt.initialChainSize = floor(NDraw*opt.initialChainSize) ;
elseif opt.initialChainSize >= NDraw
    opt.initialChainSize = NDraw ;
    utils.warning('poster:regen',...
        'Initial chain size is larger than the number of requested draws.') ;
    opt.initialChainSize = min(NDraw,initialChainSize) ;
end

% Number of estimated parameters.
nPar = length(This.paramList);

% Generate initial chain for constructing reentry distribution and special K
fprintf(1,'Generating initial chain...\n') ;
[initSample,lp_initSample,initAccRatio,initSgm,initFinalCov] ...
    = arwm(This,opt.initialChainSize,'lastAdapt=',1,'progress=',true) ; %#ok<*ASGLU,*NASGU>
initStd = chol(cov(initSample')) ;
initMean = mean(initSample,2) ;

% Construct reentry distribution
reentryDist = logdist.normal(initMean,initStd) ;
reentrySample = reentryDist([],'draw',opt.initialChainSize) ;

% Target distribution
targetDist = @(x) mylogpost(This,x,s) ;

% Construct proposal distribution
propNew = @(x) rwrand(x,chol(initFinalCov)) ;

% This special constant indirectly controls expected tour length:
% higher K means shorter expected tour length
%
% In the limit as K becomes large the algorithm reduces to a
% rejection sampling method, and as K becomes small the algorithm
% reduces to pure random walk method.
K = mean(exp(lp_initSample)) / mean(exp(reentryDist(reentrySample))) ;
ln_K = log(K)/10 ;

%--------------------------------------------------------------------------
% Main loop
NRegen = NDraw - opt.initialChainSize ;
n = 1 ;
s = 0 ;
len = zeros(1,NRegen) ;
regenSample = NaN(nPar,NRegen) ;
lp_regenSample = NaN(nPar,1) ;
Yt = [] ;
lp_Yt = [] ;
fprintf(1,'Tour     Draw     Avg Tour Length\n') ;
while n < NRegen
    [tourSample,tourLogpost,tourLength,Yt,lp_Yt] ...
        = doTour(NRegen,ln_K,nPar,reentryDist,targetDist,propNew,Yt,lp_Yt) ;
    regenSample(:,n:n+tourLength-1) = tourSample ;
    lp_regenSample(n:n+tourLength-1) = tourLogpost ;
    s = s+1 ;
    len(s) = tourLength ;
    n = n + tourLength ;
    fprintf(1,'%4.f     %4.f     %6.f\n',s,n,mean(len(1:s))) ;
end
len = len(1:s) ;

sample = [initSample regenSample] ;
lp_Sample = [lp_initSample lp_regenSample] ;


    function newTheta = rwrand(theta, sig)
        u = randn(nPar,1) ;
        newTheta = theta + sig*u ;
    end

end

function [tourSample,tourLogpost,len,Yt,lp_Yt] ...
    = doTour(N,ln_K,nPar,reentryDist,targetDist,propNew,Yt,lp_Yt)
% N: some number larger than almost any tour, affects performance but
% not results

% the Alpha state is something out of this world (or at least out of the
% support of Theta...)
alphaState = NaN(nPar,1) ;
isAlphaState = @(x) all(isnan(x),1) ; % tests which column vectors are NaN

if isempty(Yt)
    Yt = alphaState ; %start in Alpha state
    lp_Yt = NaN ;
end
tourSample = NaN(nPar,N) ;
tourLogpost = NaN(1,N) ;
len = 0 ;
accW = false ;
while ~accW
    lp_V = NaN ;
    lp_Z = NaN ;
    lp_W = NaN ;
    alpha_W = NaN ;
    accZ = false ;
    accV = false ;
    
    if isAlphaState( Yt )
        V = alphaState ;
    else
        Z = propNew( Yt ) ;
        lp_Z = targetDist( Z ) ;
        if log(rand) < min([0, lp_Z - lp_Yt])
            accZ = true ;
            [V, lp_V] = deal(Z, lp_Z) ;
        else
            [V, lp_V] = deal(Yt, lp_Yt) ;
        end
    end
    if isAlphaState( V )
        W = reentryDist([],'draw') ;
        lp_W = targetDist( W ) ;
        lq_W = reentryDist( W ) + ln_K ;
        if log(rand) < min([0, lp_W - lq_W]) ;
            accW = true ;
            [Yt, lp_Yt] = deal(W, lp_W) ;
        else
            [Yt, lp_Yt] = deal(alphaState, NaN) ;
        end
    else
        lq_V = reentryDist( V ) + ln_K ;
        if log(rand) < min([0, lq_V - lp_V])
            accV = true ;
            [Yt, lp_Yt] = deal(V, lp_V) ;
        else
            [Yt, lp_Yt] = deal(alphaState, NaN) ;
        end
    end
    
    if isAlphaState( Yt )
        % don't store alpha states
    else
        len = len + 1 ;
        tourSample(:,len) = Yt ;
        tourLogpost(len) = lp_Yt ;
    end
end

tourSample = tourSample(:,1:len) ;
tourLogpost = tourLogpost(1:len) ;

end