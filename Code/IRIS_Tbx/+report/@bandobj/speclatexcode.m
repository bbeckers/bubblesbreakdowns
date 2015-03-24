function [C,Temps] = speclatexcode(This)
% speclatexcode  [Not a public function] LaTeX code for report/band data.
%
% Backend IRIS function.
% No help provided.

% -IRIS Toolbox.
% -Copyright (c) 2007-2013 IRIS Solutions Team.

%--------------------------------------------------------------------------

Temps = {};

par = This.parent;
time = par.options.range;
colStruct = par.options.colstruct;

[cData,time] = getdata(This,This.data,time,colStruct);
lData = getdata(This,This.low,time,colStruct);
hData = getdata(This,This.high,time,colStruct);
cData = cData(:,:);
lData = lData(:,:);
hData = hData(:,:);

text = This.caption;
br = sprintf('\n');
C = '';

nx = max([size(cData,2),size(lData,2),size(hData,2)]);
for iRow = 1 : nx
    iCData = cData(:,min(iRow,end));
    iLData = lData(:,min(iRow,end));
    iHData = hData(:,min(iRow,end));
    if This.options.relative
        iLData = iCData + iLData;
        iHData = iCData + iHData;
    end
    iData = [iCData,iLData,iHData];
    if iRow <= numel(This.options.marks)
        mark = This.options.marks{iRow};
    else
        mark = '';
    end
    C = [C,br, ...
        latexonerow(This,iRow,time,iData,mark,text)]; ...
        %#ok<AGROW>
end

end