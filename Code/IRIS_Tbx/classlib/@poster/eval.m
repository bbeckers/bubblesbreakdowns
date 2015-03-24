function [Obj,L,PP,SP] = eval(This,varargin)
% eval  Evaluate posterior density at specified points.
%
% Syntax
% =======
%
%     [X,L,PP,SrfP,FrfP] = eval(Pos)
%     [X,L,PP,SrfP,FrfP] = eval(Pos,P)
%
% Input arguments
% ================
%
% * `Pos` [ poster ] - Posterior object returned by the
% [`model/estimate`](model/estimate) function.
%
% * `P` [ struct ] - Struct with parameter values at which the posterior
% density will be evaluated; if `P` is not specified, the posterior density
% at the point of the estimated mode is returned.
%
% Output arguments
% =================
%
% * `X` [ numeric ] - The value of log posterior density evaluated at `P`;
% N.B. the returned value is log posterior, and not minus log posterior.
%
% * `L` [ numeric ] - Contribution of data likelihood to log posterior.
%
% * `PP` [ numeric ] - Contribution of parameter priors to log posterior.
%
% * `SrfP` [ numeric ] - Contribution of shock response function priors to
% log posterior.
%
% * `FrfP` [ numeric ] - Contribution of frequency response function priors
% to log posterior.
%
% Description
% ============
%
% The total log posterior consists, in general, of the four contributions
% listed above:
%
%     X = L + PP + SrfP + FrfP.
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

if isempty(varargin)
    p = This.initParam;
elseif length(varargin) == 1
    p = varargin{1};
else
    p = varargin;
end

%--------------------------------------------------------------------------

if nargin == 1 && nargout <= 1
    % Return log posterior at optimum.
    Obj = This.initLogPost;
    return
end

s = mylogpoststruct(This);

% Evaluate log poeterior at specified parameter sets. If
% it's multiple parameter sets, pass them in as a cell, not
% as multiple input arguments.
if isstruct(p)
    p0 = p;
    nPar = length(This.paramList);
    p = nan(1,nPar);
    for i = 1 : nPar
        p(i) = p0.(This.paramList{i});
    end
end

if ~iscell(p)
    p = {p};
end
np = numel(p);

% Minus log posterior.
Obj = nan(size(p));
% Minus log likelihood.
L = nan(size(p));
% Minus log parameter priors.
PP = nan(size(p));
% Minus log system priors.
SP = nan(size(p));

parfor i = 1 : np
    theta = p{i}(:);
    [Obj(i),L(i),PP(i),SP(i)] = mylogpost(This,theta,s);
end

end