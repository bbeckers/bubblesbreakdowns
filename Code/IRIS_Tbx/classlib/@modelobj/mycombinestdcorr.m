function Stdcorr = mycombinestdcorr(ThisStdcorr,UsrStdcorr,NPer)
% mycombinestdcorr  [Not a public function] Combine model stdcorr vector with user-supplied time-varying stdcorr.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

usrStdcorrInx = ~isnan(UsrStdcorr);
ThisStdcorr = ThisStdcorr(:);
if any(usrStdcorrInx(:))
    lastUser = max(1,size(UsrStdcorr,2));
    Stdcorr = ThisStdcorr(:,ones(1,lastUser));
    Stdcorr(usrStdcorrInx) = UsrStdcorr(usrStdcorrInx);
    % Add the model stdcorrs if the last user-supplied data point is before
    % the end of the sample.
    if size(Stdcorr,2) < NPer
        Stdcorr(:,end+1) = ThisStdcorr;
    end
else
    Stdcorr = ThisStdcorr;
end

end