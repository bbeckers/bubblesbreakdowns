function irisset(varargin)
% irisset  Change configurable IRIS options.
%
% Syntax
% =======
%
%     irisset(Option,Value)
%     irisset(Option,Value,Option,Value,...)
%
% Input arguments
% ================
%
% * `Option` [ char ] - Name of the IRIS configuration option that will be
% modified.
%
% * `Value` [ ... ] - New value that will be assigned to the option.
%
% Modifiable IRIS configuration options
% ======================================
%
% Dates and formats
% -------------------
%
% * `'dateFormat='` [ char | *'YPF'* ] - Date format used to display dates
% in the command window, CSV databases, and reports. Note that the default
% date format for graphs is controlled by the `'plotdateformat='` option.
% The default 'YFP' means that the year, frequency letter, and period is
% displayed. See also help on [`dat2str`](dates/dat2str) for more date
% formatting details. The `'dateformat='` option is also found in many IRIS
% functions whenever it is relevant, and can be used to overwrite the
% `'irisset='` settings.
%
% * `'freqLetters='` [ char | *'YHQBM'* ] - Five letters used to represent
% the five possible frequencies of the IRIS dates: yearly, half-yearly,
% quarterly, bi-monthly, and monthly, such as the 'Q' in '2010Q1' denoting.
%
% * `'months='` [ cellstr | *{'January',...,'December'}* ] - Twelve strings
% representing the names of the twelve months. This option can be used
% whenever you want to replace the default English names with your local
% language. .
%
% * `'plotDateFormat='` [ char | *'Y:P'* ] - Date format used to display dates in
% graphs including graphs in reports. The default is 'Y:P', which is
% different from the `'dateformat='` option.
%
% * `'tseriesFormat='` [ char | *empty* ] - Format string for displaying time
% series data on the screen. See help on the Matlab `sprintf` function for
% how to set up format strings. If empty the default format of the
% `num2str` function is used.
%
% * `tseriesMaxWSpace='` [ numeric | *`5`* ] - Maximum number of white
% spaces printed between individual columns of a multivariate tseries
% object on the screen.
%
% * `'standinMonth='` [ *'first'* | 'last' | numeric ] - This option
% specifies which month will be used to represent lower-frequency periods
% (such as a quarters) when a month-displaying format is used in
% `'dateformat='`.
%
% External tools used by IRIS
% -----------------------------
%
% * `'pdflatexPath='` [ char ] - Location of the `pdflatex.exe` program. This
% program is called to compile report and publish m-files. By
% default, IRIS attempts to locate `pdflatex.exe` by running TeX's
% `kpsewhich`, and `which` on Unix platforms.
%
% * `'epstopdfPath='` [ char ] - Location of the `epstopdf.exe` program.
% This program is called to convert EPS graphics files to PDFs in
% reports.
%
% Other options
% ---------------
%
% * `'extensions='` [ cellstr | *{'model','s','q'}* ] - List of extensions
% that are automatically associated with the Matlab editor. The advantage
% of such associations is that the IRIS [model files](modellang/Contents),
% [steady-state files](sstatelang/Contents), and
% [quick-report](qreport/Contents) files get syntax-highlighted.
%
% Description
% ============
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

irisconfigmaster('set',varargin{:});

end