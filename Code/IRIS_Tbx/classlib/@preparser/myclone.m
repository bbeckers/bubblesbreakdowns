function C = myclone(C,Clone)
% myclone  [Not a public function] Clone a preparsed code.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

if ~preparser.mychkclonestring(Clone)
    utils.error('preparser', ...
        'Invalid clone string: ''%s''.', ...
        Clone);
end

pattern = '(?<!!)\<([A-Za-z]\w*)\>(?!\()';
replace = '${strrep(Clone,''?'',$0)}';
C = regexprep(C,pattern,replace);

end