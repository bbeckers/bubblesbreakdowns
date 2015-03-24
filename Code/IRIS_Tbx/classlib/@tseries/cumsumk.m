function X = cumsumk(X,K,RHO,RANGE,varargin)
% cumsumk  Cumulative sum with a k-period leap.
%
% Syntax
% =======
%
%     Y = cumsumk(X,K,RHO,RANGE)
%     Y = cumsumk(X,K,RHO)
%     Y = cumsumk(X,K)
%     Y = cumsumk(X)
%
% Input arguments
% ================
%
% * `X` [ tseries ] - Input data.
%
% * `K` [ numeric ] - Number of periods that will be leapt the cumulative
% sum will be taken; if not specified, `K` is chosen to match the frequency
% of the input data (e.g. `K = -4` for quarterly data), or `K = -1` for
% indeterminate frequency.
%
% * `RHO` [ numeric ] - Autoregressive coefficient; if not specified, `RHO
% = 1`.
%
% * `RANGE` [ numeric ] - Range on which the cumulative sum will be
% computed and the output series returned.
%
% Output arguments
% =================
%
% * `Y` [ tseries ] - Output data constructed as described below.
%
% Options
% ========
%
% * `'log='` [ `true` | *`false`* ] - Logarithmise the input data before,
% and de-logarithmise the output data back after, running `x12`.
%
% Description
% ============
%
% If `K < 0`, the first `K` observations in the output series `Y` are
% copied from `X`, and the new observations are given recursively by
%
%     Y{t} = RHO*Y{t-K} + X{t}.
%
% If `K > 0`, the last `K` observations in the output series `Y` are
% copied from `X`, and the new observations are given recursively by
%
%     Y{t} = RHO*Y{t+K} + X{t},
%
% going backwards in time.
%
% If `K == 0`, the input data are returned.
%
% Example
% ========
%
% Construct random data with seasonal pattern, and run X12 to seasonally
% adjust these series.
%
%     x = tseries(qq(1990,1):qq(2020,4),@randn);
%     x1 = cumsumk(x,-4,1);
%     x2 = cumsumk(x,-4,0.7);
%     x1sa = x12(x1);
%     x2sa = x12(x2);
%
% The new series `x1` will be a unit-root process while `x2` will be
% stationary. Note that the command on the second line could be replaced
% with `x1 = cumsumk(x)`.
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

try
    K;
catch
    K = -max(1,datfreq(X.start));
end         

try
    RHO;
catch
    RHO = 1;
end

try
    RANGE; %#ok<*VUNUS>
catch %#ok<*CTCH>
    RANGE = Inf;
end

opt = passvalopt('tseries.cumsumk',varargin{:});

if K == 0
    return
end

%**************************************************************************

datasize = size(X.data);
X.data = X.data(:,:);
[data,range] = rangedata(X,RANGE);

if opt.log
    data = log(data);
end

data = tseries.mycumsumk(data,K,RHO);

if opt.log
    data = exp(data);
end

X.start = range(1);
if length(datasize) == 2
    X.data = data;
else
    X.data = reshape(data,[size(data,1),datasize(2:end)]);
end

end