% Navigate to Data folder
cd('..\Daten\Realtime, monatlich als csv')
addpath('..\..\Code\Functions')

clear
clc

%% Find and select data files
list = dir;
list2 = cell(length(list)-2,1);
for i=3:length(list)
    list2{i-2} = list(i,1).name;
end

del = [];
for i=1:length(list2)
    if ~strcmp(list2{i}(end-2:end),'csv')
        del = [del;i];
    end
end
list2(del) = [];

list = list2;
nfiles = length(list);
varnames = cell(nfiles,1);
for i=1:nfiles
    varnames{i} = list{i}(1:end-4);
end
clear list2 del i dir

%% Load data and perform tests
% Maximum lag length
pmax = 12;

% Outputs
hADF_notrend = zeros(nfiles,1);
hADF_trend = zeros(nfiles,1);
pValADF_notrend1 = zeros(nfiles,1);
pValADF_trend1 = zeros(nfiles,1);
pValADF_notrend2 = zeros(nfiles,1);
pValADF_trend2 = zeros(nfiles,1);
hKPSS_notrend = zeros(nfiles,1);
hKPSS_trend = zeros(nfiles,1);
pValKPSS_notrend1 = zeros(nfiles,1);
pValKPSS_trend1 = zeros(nfiles,1);
pValKPSS_notrend2 = zeros(nfiles,1);
pValKPSS_trend2 = zeros(nfiles,1);
pstarnotrend1 = zeros(nfiles,1);
pstartrend1 = zeros(nfiles,1);
pstarnotrend2 = zeros(nfiles,1);
pstartrend2 = zeros(nfiles,1);
trend1 = zeros(nfiles,1);
trendpVal1 = zeros(nfiles,1);
trend2 = zeros(nfiles,1);
trendpVal2 = zeros(nfiles,1);
hADF = zeros(pmax+1,2,nfiles);
pValADF = zeros(pmax+1,2,nfiles);
hKPSS = zeros(pmax+1,2,nfiles);
pValKPSS = zeros(pmax+1,2,nfiles);

unitroot1 = zeros(nfiles,1);
unitroot2 = zeros(nfiles,1);

integlevel = zeros(nfiles,4);

for i=1:nfiles
    [data,vintages] = xlsread(list{i});
    vintages = vintages(1,:);
    y = data(:,end);
    y = log(y(~isnan(y)));
    %% Stationarity Tests
    [hADF_notrend(i),hADF_trend(i),pValADF_notrend1(i),pValADF_trend1(i),...
    hKPSS_notrend(i),hKPSS_trend(i),pValKPSS_notrend1(i),pValKPSS_trend1(i),...
    pstarnotrend1(i),pstartrend1(i),trend1(i),trendpVal1(i),...
    hADF(:,:,i),pValADF(:,:,i),hKPSS(:,:,i),pValKPSS(:,:,i)]...
    = stationarity(y,pmax);
    % Integration level
    if hADF_notrend(i)==0;
        integlevel(i,1) = integlevel(i,1)+1;
    elseif hADF_trend(i)==0;
        integlevel(i,2) = integlevel(i,2)+1;
    elseif hKPSS_notrend(i)==1;
        integlevel(i,3) = integlevel(i,3)+1;
    elseif hKPSS_trend(i)==1;
        integlevel(i,4) = integlevel(i,4)+1;
    end
    unitroot1(i) = sum([1-hADF_notrend;1-hADF_trend;hKPSS_notrend;hKPSS_trend])>0;
    if unitroot1(i)==1
        dy = y(2:end)-y(1:end-1);
        [hADF_notrend(i),hADF_trend(i),pValADF_notrend2(i),pValADF_trend2(i),...
        hKPSS_notrend(i),hKPSS_trend(i),pValKPSS_notrend2(i),pValKPSS_trend2(i),...
        pstarnotrend2(i),pstartrend2(i),trend2(i),trendpVal2(i),...
        hADF(:,:,i),pValADF(:,:,i),hKPSS(:,:,i),pValKPSS(:,:,i)]...
        = stationarity(dy,pmax);
        % Integration level
        if hADF_notrend(i)==0;
            integlevel(i,1) = integlevel(i,1)+1;
        elseif hADF_trend(i)==0;
            integlevel(i,2) = integlevel(i,2)+1;
        elseif hKPSS_notrend(i)==1;
            integlevel(i,3) = integlevel(i,3)+1;
        elseif hKPSS_trend(i)==1;
            integlevel(i,4) = integlevel(i,4)+1;
        end
        unitroot2(i) = sum([1-hADF_notrend;1-hADF_trend;hKPSS_notrend;hKPSS_trend])>0;
        if unitroot2(i)==1
            dy2 = dy(2:end)-dy(1:end-1);
        end
    end
end

cd('..\..\Code')

xlswrite('Stationarity_Analysis.xlsx',varnames,'Integration','A3')
xlswrite('Stationarity_Analysis.xlsx',integlevel,'Integration','B3')

xlswrite('Stationarity_Analysis.xlsx',varnames,'Trend','A2')
xlswrite('Stationarity_Analysis.xlsx',trend1,'Trend','B2')
xlswrite('Stationarity_Analysis.xlsx',trendpVal1,'Trend','C2')

xlswrite('Stationarity_Analysis.xlsx',varnames,'Lag length','A2')
xlswrite('Stationarity_Analysis.xlsx',pstarnotrend1,'Lag length','B2')
xlswrite('Stationarity_Analysis.xlsx',pstarnotrend2,'Lag length','C2')