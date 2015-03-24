function Flag = rngcmp(V1,V2)
% rngcmp  True if two VAR objects have been estimated using the same dates.
%
% Syntax
% -------
%
%     Flag = rngcmp(V1,V2)
%
% Input arguments
% ================
%
% * `V1`, `V2` [ VAR ] - Two estimated VAR objects.
%
% Output arguments
% =================
% 
% * `Flag` [ `true` | `false` ] - True if the two VAR objects, `V1` and
% `V2`, have been estimated using observations at the same dates.
%
% Description
% ============
%
% Example
% ========
%

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

pp = inputParser();
pp.addRequired('V1',@(x) isa(x,'VAR'));
pp.addRequired('V2',@(x) isa(x,'VAR'));
pp.parse(V1,V2);

%--------------------------------------------------------------------------

nAlt1 = size(V1.A,3);
nAlt2 = size(V2.A,3);
nAlt = max(nAlt1,nAlt2);

Flag = false(1,nAlt);
for iAlt = 1 : nAlt
    fitted1 = V1.fitted(:,:,min(iAlt,end));
    fitted2 = V2.fitted(:,:,min(iAlt,end));
    range1 = V1.range(fitted1);
    range2 = V2.range(fitted2);
    Flag(iAlt) = length(range1) == length(range2) ...
        && all(datcmp(range1,range2));
end

end