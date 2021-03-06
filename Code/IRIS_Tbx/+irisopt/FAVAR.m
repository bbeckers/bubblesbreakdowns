function def = FAVAR()
% FAVAR   [Not a public function] Default options for FAVAR class functions.

%**************************************************************************

def = struct();

def.estimate = { ...
    'cross',true, ...
    @(x) islogicalscalar(x) || (isnumericscalar(x) && x >=0 && x <= 1), ...
    'method','auto',@(x) isequal(x,'auto') || isequal(x,1) || isequal(x,2), ...
    'order',1,@(x) isnumericscalar(x), ...
    'output','auto',@(x) ischar(x) && any(strcmpi(x,{'auto','dbase','tseries','array'})), ...
    'rank',Inf,@(x) isnumericscalar(x), ...
    'tolerance','auto',@(x) strcmpi(x,'auto') || isnumericscalar(x), ...
    'ynames,yname',@(n) ['y',sprintf('%g',n)],@(x) iscellstr(x) || isfunc(x), ...
    };

def.filter = { ...
    'cross',true, ...
    @(x) islogicalscalar(x) || (isnumericscalar(x) && x >=0 && x <= 1), ...
    'invfunc','auto',@(x) isequal(x,'auto') || isfunc(x), ...
    'meanonly',false,@islogicalscalar, ...
    'output','auto',@(x) ischar(x) && any(strcmpi(x,{'auto','dbase','tseries','array'})), ...
    'persist',false,@islogicalscalar, ...
    'tolerance',0,@(x) isnumericscalar(x), ...
    };

def.forecast = { ...
    'cross',true, ...
    @(x) islogicalscalar(x) || (isnumericscalar(x) && x >=0 && x <= 1), ...
    'invfunc','auto',@(x) isequal(x,'auto') || isfunc(x), ...
    'meanonly',false,@islogicalscalar, ...
    'output','auto',@(x) ischar(x) && any(strcmpi(x,{'auto','dbase','tseries','array'})), ...
    'persist',false,@islogicalscalar, ...
    'tolerance',0,@(x) isnumericscalar(x), ...
    };

end