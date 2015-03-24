function flag = issolved(m)
% issolved  True if a model solution exists.
%
% Syntax
% =======
%
%     flag = issolved(m)
%
% Input arguments
% ================
%
% * `m` [ model ] - Model object.
%
% Output arguments
% =================
%
% * `flag` [ `true` | `false` ] - True for each parameterisation for which a
% stable unique solution has been found and exists currently in the model
% object.
%
% Description
% ============
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%**************************************************************************

[ans,flag] = isnan(m,'solution');
flag = ~flag;

end
