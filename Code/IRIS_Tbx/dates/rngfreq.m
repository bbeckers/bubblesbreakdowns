function Freq = rngfreq(Range)
% rngfreq  Date frequency of a date range.
%
% Backed IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

if iscell(Range)
    Freq = nan(size(Range));
    for i = 1 : numel(Range)
        Freq(i) = rngfreq(Range{i});
    end
    return
end

%--------------------------------------------------------------------------

Range = Range(:).';

if any(isnan(Range))
    Freq = NaN;
    return
end

if all(isinf(Range))
    Freq = Inf;
    return
end

Range = Range(~isinf(Range));

freq = datfreq(Range);
if all(freq == freq(1))
    Freq = freq(1);
else
    Freq = NaN;
end

end