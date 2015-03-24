function Def = tseries()
% tseries  [Not a public function] Default options for tseries class functions.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

Def = struct();

convert = { ...
    'function',[],@(x) isempty(x) || isfunc(x) || ischar(x), ...
    'missing',NaN,@(x) (ischar(x) && any(strcmpi(x,{'last'}))) || isnumericscalar(x), ...
    };

config = irisget();

dates = { ...
    'dateformat','config',config.validate.dateformat, ...
    'freqletters,freqletter','config',config.validate.freqletters, ...
    'months,month','config',config.validate.months, ...
    'standinmonth','config',config.validate.standinmonth, ...
    'datetick,dateticks',Inf,@(x) isnumeric(x) ...
    || any(strcmpi(x,{'yearstart','yearend','yearly'})) ...
    || isfunc(x), ...
    };

Def.acf = { ...
    'demean',true,@islogicalscalar, ...
    'order',0,@isnumericscalar, ...
    'smallsample',true,@islogicalscalar, ...
    };

Def.bpass = { ...
    'addtrend',true,@islogicalscalar, ...
    'detrend',true,@islogicalscalar, ...
    'log',false,@islogicalscalar, ...
    'method','cf',@(x) ischar(x) && any(strcmpi(x,{'cf','hwfsf'})), ...
    'ttrend',[],@(x) isempty(x) || islogicalscalar(x), ...
    'unitroot',true,@islogicalscalar, ...
    'window','hamming',@(x) ischar(x) && any(strcmpi(x,{'hamming','hanning','none'})), ...
    };

Def.chowlin = { ...
    'constant',true,@islogicalscalar, ...
    'log',false,@islogicalscalar, ...
    'ngrid',200,@(x) isnumericscalar(x) && x >= 1, ...
    'rho','estimate', ...
    @(x) any(strcmpi(x,{'auto','estimate','negative','positive'})) ...
    || (isnumericscalar(x) && x > -1 && x < 1), ...
    'timetrend',false,@islogicalscalar, ...
    };

Def.convertdaily = [convert,{ ...
    'ignorenan',true,@islogical, ...
    'method',@mean,@(x) isfunc(x) || isequal(x,'first') || isequal(x,'last'), ...
    'select',Inf,@(x) isnumeric(x), ...
    }];

Def.convertaggregate = [convert,{ ...
    'ignorenan',false,@islogical, ...
    'method',@mean,@(x) isfunc(x) || isequal(x,'first') || isequal(x,'last'), ...
    'select',Inf,@(x) isnumeric(x), ...
    }];

Def.convertinterp = [convert,{ ...
    'ignorenan',false,@islogical, ...
    'method','cubic',@(x) ischar(x), ...
    'position','centre',@(x) ischar(x) && any(strncmpi(x,{'c','s','e'},1)), ...
    }];

Def.cumsumk = { ...
    'log',false,@islogicalscalar, ...
    };

Def.errorbar = { ...
    'relative',true,@islogicalscalar, ...
    };

Def.expsmooth = { ...
    'init',NaN,@isnumeric, ...
    'log',false,@islogicalscalar, ...
};

Def.fft = { ...
    'full',false,@islogicalscalar, ...
    };

Def.filter = { ...
    'change,growth',[],@(x) isempty(x) || istseries(x), ...
    'gamma',1,@(x) isa(x,'tseries') || (isnumericscalar(x) && x > 0), ...
    'cutoff',[],@(x) isempty(x) || (isnumeric(x) && all(x(:) > 0)), ...
    'cutoffyear',[],@(x) isempty(x) || (isnumeric(x) && all(x(:) > 0)), ...
    'drift',0,@(x) isnumeric(x) && length(x) == 1, ...
    'gap',[],@(x) isempty(x) || istseries(x), ...
    'infoset',2,@(x) isequal(x,1) || isequal(x,2), ...
    'lambda',[],@(x) isempty(x) || (isnumeric(x) && all(x(:) > 0)) || (ischar(x) && strcmpi(x,'auto')), ...
    'level',[],@(x) isempty(x) || istseries(x), ...
    'log',false,@islogical, ...
    'swap',false,@islogical, ...
    'forecast',[],@(x) isnumeric(x) && length(x) <= 1, ...
    };

Def.myplot = { ...
    dates{:}, ...
    'dateposition','c',@(x) ischar(x) && ~isempty(x) && any(x(1) == 'sec'), ...
    'function',[],@(x) isempty(x) || isfunc(x), ...
    'tight',false,@islogicalscalar, ...
    'xlimmargin','auto',@(x) islogicalscalar(x) || (ischar(x) && strcmpi(x,'auto')), ...
    }; %#ok<CCAT>

Def.interp = { ...
    'method','cubic',@ischar, ...
    };

Def.moving = { ...
    'window',Inf,@isnumeric, ...
    'function',@mean,@isfunc, ...
    };

Def.barcon = {
    'barwidth',0.8,@isnumericscalar, ...
    'colormap',[],@isnumeric, ...
    'evenlyspread',true,@islogicalscalar, ...
    'ordering','preserve',@(x) isanychari(x,{'descend','ascend','preserve'}) ...
    || isnumeric(x), ...
    };

Def.normalise = { ...
    'mode','mult',@(x) any(strncmpi(x,{'add','mult'},3)), ...
    };

Def.pct = { ...
    'outputfreq,freq',Inf,@(x) isnumericscalar(x) && any(x == [Inf,1,2,4,6,12]), ...
    };

Def.plotcmp = { ...
    'compare',[-1;1],@isnumeric, ...
    'cmpcolor,diffcolor',[1,0.75,0.75],@(x) isnumeric(x) && length(x) == 3 && all(x >= 0) && all(x <= 1), ...
    'rhsplotfunc',[],@(x) isempty(x) || isanychari(char(x),{'bar','area'}), ...
    'cmpplotfunc,diffplotfunc',@bar,@(x) isanychari(char(x),{'bar','area'}), ...
    };

Def.plotpred = { ...
    'connect',true,@islogicalscalar, ...
    'firstline',{},@(x) iscell(x) && iscellstr(x(1:2:end)), ...
    'predlines',{},@(x) iscell(x) && iscellstr(x(1:2:end)), ...
    'firstmarker,firstmarkers,startpoint,startpoints','.',@(x) ischar(x) ...
    && any(strcmpi(x,{'none','+','o','*','.','x','s','d','^','v','>','<','p','h'})), ...
    'shownanlines',true,@islogicalscalar, ...
};

Def.plotyy = { ...
    Def.myplot{:}, ...
    'coincident',false,@islogicalscalar, ...
    'highlight',[],@isnumeric, ...
    'lhsplotfunc',@plot,@(x) ischar(x) || isfunc(x), ...
    'rhsplotfunc',@plot,@(x) ischar(x) || isfunc(x), ...
    'lhstight',false,@islogicalscalar, ...
    'rhstight',false,@islogicalscalar, ...
    };

Def.regress = {...
    'constant',false,@islogical, ...
    'weighting',[],@(x) (isnumeric(x) && isempty(x)) || istseries(x),...
    };

Def.spy = { ...
    Def.myplot{:}, ...
    'names,name',{},@iscellstr, ...
    'test',@isfinite,@(x) isfunc(x), ...    
    };

Def.trend = { ...
    'break',[],@isnumeric, ...
    'connect',true,@islogicalscalar, ...
    'diff',false,@islogicalscalar, ...
    'log',false,@islogicalscalar, ...
    'season',false,@(x) isempty(x) || islogicalscalar(x) || isnumericscalar(x), ...
    };

Def.windex = { ...
    'log',false,@islogical, ...
    'method','simple',@(x) ischar(x) && any(strcmpi(x,{'simple','divisia'})), ...
    };

Def.x12 = { ...
    'arima',[], @(x) isempty(x) || islogicalscalar(x) ...
    || (isnumeric(x) && length(x) <= 2 && all(x >= 0) && all(x == round(x))), ...
    'backcast,backcasts',0,@(x) isnumericscalar(x), ...
    'cleanup,deletetempfiles,deletetempfile,deletex12file,deletex12file,delete',true,@islogicalscalar, ...
    'dummy',[],@(x) isempty(x) || isa(x,'tseries'), ...
    'dummytype','holiday',@(x) ischar(x) && any(strcmpi(x,{'holiday','td','ao'})), ...
    'display',false,@islogicalscalar, ...
    'forecast,forecasts',0,@(x) isnumericscalar(x), ...
    'log',false,@islogicalscalar, ...
    'maxiter',1500,@(x) isnumericscalar(x) && x > 0 && x == round(x), ...
    'maxorder',[2,1],@(x) isnumeric(x) && length(x) == 2 ...
    && any(x(1) == [1,2,3,4]) && any(x(2) == [1,2]), ...
    'missing',false,@islogicalscalar, ...
    'mode','auto',@(x) (isnumeric(x) && any(x == -1 : 3)) || any(strcmp(x,{'add','a','mult','m','auto','sign','pseudo','pseudoadd','p','log','logadd','l'})), ...
    'output','d11',@(x) ischar(x) || iscellstr(x), ...
    'saveas','',@ischar, ...
    'specfile','default',@(x) ischar(x) || isinf(x), ...
    'tdays,tday',false,@islogicalscalar, ...
    'tolerance',1e-5,@(x) isnumericscalar(x) && x > 0, ...
    };

end