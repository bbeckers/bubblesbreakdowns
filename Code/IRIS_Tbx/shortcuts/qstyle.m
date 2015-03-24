function varargout = qstyle(varargin)
% qplot  Shortcut for qreport.qstyle.
%
% See help on [`qreport.qstyle`](qreport/qstyle).

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

[varargout{1:nargout}] = qreport.qstyle(varargin{:});

end