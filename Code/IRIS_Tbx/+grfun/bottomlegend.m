function Le = bottomlegend(varargin)
% bottomlegend  Horizontal graph legend displayed at the bottom of the figure window.
%
% Syntax
% =======
%
%     Le = grfun.bottomlegend(Entry,Entry,...)
%
% Input arguments
% ================
%
% * `Entry` [ char | cellstr ] - Legend entries; same as in the standard
% `legend` function.
%
% Output arguments
% =================
%
% * `AX` [ numeric ] - Handle to the legend axes object created.
%
% Description
% ============
%
% Example
% ========

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

Le = legend(varargin{:});
if ~isempty(Le)
    set(Le,'orientation','horizontal');
    Le = grfun.movetosubplot(Le,'bottom');
end

end