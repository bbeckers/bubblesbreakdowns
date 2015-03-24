function M = per2month(Per,Freq,StandinMonth)
% per2month  [Not a public function] Return month to represent a given period.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

if ischar(StandinMonth)
    switch StandinMonth
        case {'first','start'}
            StandinMonth = 1;
        case {'last','end'}
            StandinMonth = 12/Freq;
        otherwise
            StandinMonth = 1;
    end
end

M = (Per-1).*12./Freq + StandinMonth;

end