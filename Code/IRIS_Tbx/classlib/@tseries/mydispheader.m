function mydispheader(This)
% mydispheader  [Not a public function] Display header for tseries objects.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

tmpSize = size(This.data);
nper = tmpSize(1);
fprintf('\t');
if isempty(This.data)
   fprintf('empty ');
end
fprintf('tseries object: %g%s\n',nper,sprintf('-by-%g',tmpSize(2:end)));
strfun.loosespace();

end