function [This,Data,B,Count] = SVAR(V,Data,varargin)
% SVAR  [Not a public function] Identify SVAR from a reduced-form VAR.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

if nargin == 1
    Data = [];
end

% Parse required input arguments.
pp = inputParser();
pp.addRequired('V',@(x) isa(x,'VAR'));
pp.addRequired('Data',@(x) isempty(x) || isnumeric(x) || istseries(x) ...
    || isstruct(x));
pp.parse(V,Data);

opt = passvalopt('SVAR.SVAR',varargin{1:end});

%--------------------------------------------------------------------------

ny = size(V.A,1);
nAlt = size(V.A,3);

% Create an empty SVAR object.
This = SVAR();
This.B = nan(ny,ny,nAlt);
This.std = nan(1,nAlt);

% Populate the superclass VAR properties.
list = utils.ndprop('VAR');
nList = length(list);
for i = 1 : nList
   This.(list{i}) = V.(list{i});
end

% Identify the B matrix.
[This,Data,B,Count] = myidentify(This,Data,opt);

end