function Flag = isempty(This,varargin)
% isempty  True if the system priors object is empty.
%
% Syntax
% =======
%
%     Flag = isempty(S)
%
% Input arguments
% ================
%
% * `S` [ systempriors ] - System priors,
% [`systempriors`](systempriors/Contents), object.
%
% Output arugments
% =================
%
% * `Flag` [ true | false ] - True if the system priors object, `S`, is
% empty, false otherwise.
%
% Description
% ============
%
% Example
% ========

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

if isempty(varargin)
    Flag = isempty(This.eval);
else
    Flag = isempty(This.systemFunc.(lower(varargin{1})).page);
end

end