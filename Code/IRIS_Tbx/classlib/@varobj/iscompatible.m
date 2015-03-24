function Flag = iscompatible(V1,V2)
% iscompatible  [Not a public function] True if two varobj objects can occur together on the LHS and RHS in an assignment.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

try
    Flag = isa(V1,'varobj') && isa(V2,'varobj') ...
        && isequal(V1.Ynames,V2.Ynames) ...
        && isequal(V1.Enames,V2.Enames) ...
        && rngcmp(V1.range,V2.range) ...
        && size(V1.A,1) == size(V2.A,1) ...
        && size(V1.A,2) == size(V2.A,2) ...
        && size(V1.K,1) == size(V2.K,1) ...
        && size(V1.K,2) == size(V2.K,2);        
catch %#ok<CTCH>
    Flag = false;
end

end