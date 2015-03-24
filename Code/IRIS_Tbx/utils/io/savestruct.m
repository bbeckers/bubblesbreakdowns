function savestruct(FName,X)
% savestruct  Convert an object to a struct and save field by field.
%
% Syntax
% =======
%
%     savestruct(FName,X)
%
% Input arguments
% ================
%
% * `FName` [ char ] - File name.
%
% * `X` [ .... ] - Object to be saved.
%
% Description
% ============
%
% `savestruct` and `loadstruct` were introduced to deal with some
% inefficiencies in standard saving and loading procedures in older
% Matlabs. In current versions of Matlab, this is no longer necessary, and
% `savestruct` and `loadstruct` functions are considered obsolete.
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

% Allow both savestruct(fname,d) and savestruct(d,fname).
if (isobject(FName) || isstruct(FName) || iscell(FName)) && ischar(X)
    [FName,X] = deal(X,FName);
end

cl = class(X);
switch cl
    case {'model','VAR','container','tseries'}
        status = warning();
        warning('off','MATLAB:structOnObject');
        X = saveobj(X);
        X = struct(X);
        warning(status);
    case 'struct'
        % Do nothing.
    otherwise
        utils.error('Cannot save %s objects using SAVESTRUCT.',cl);
end

X.SAVESTRUCT_CLASS = cl; %#ok<STRNU>

% Save individual fields of underlying struct.
save(FName,'-struct','x','-mat');

end