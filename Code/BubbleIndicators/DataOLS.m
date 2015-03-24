%% Load Data
cd(strcat('..\Data\',country))
addpath(genpath('..\..\Matlab'));
% Load non-revised data
D = dbload('NonrevData.csv','nameRow=',1,'commentrow=','Descriptor');

% Load real-time data
cd('Real-time') % Choose folder
files = dir('*csv*'); % List of all csv-files
M = cell(2,1); % Load files in cells
for i=1:length(files)
    trash = dbload(files(i,1).name,'nameRow=',1);
    M = [M(1,:), fieldnames(trash)'; M(2,:), struct2cell(trash)'];
end; clear trash files i
M = M(:,2:end);
DRT = struct(M{:}); % Transform to structure
clear M
vintages = fieldnames(DRT);

cd('..\..\..\Matlab')

% Set new base years
date.base = (mm(1995,1):mm(1995,12))';
date.baseq = (qq(1995,1):qq(1995,4))';

%% Lag structure, corrections and variable declaration
% Publication lags
l = struct();
l.gdp = (1:3)';
l.defl = (1:3)';
l.cpi = 1;
l.uc = 1;
l.ip = 1;

% Number of corrections
c = struct();
c.gdp = 2;
c.defl = 2;
c.cpi = 0;
c.uc = 0;
c.ip = 2;

% Variables loaded from database
vars = {'awh', 'cpi', 'defl', 'emp', 'gdp', 'ip', 'renti', 'uc'};
% Real-time variables
varsrt = {'awh', 'defl', 'emp', 'gdp', 'renti', 'uc', 'ip'};

%% Transformations
% Historical database
C = cell(4,length(fieldnames(D)));
C(1:2,:) = [fieldnames(D)'; struct2cell(D)'];
for i=1:length(fieldnames(D))
    C{3,i} = lower(C{1,i}); % Convert series names
    % Take h-step log-differences of all time-series except interest rates
    if strcmp(C{1,i},'CPI')==1 || strcmp(C{1,i},'IP')==1
        C{4,i} = NaN(size(C{2,i},1),size(par.h,1));
        startvint = get(C{2,i},'start');
        endvint = get(C{2,i},'end');
        for h=1:size(par.h,1)
            if strcmp(C{1,i},'CPI')
                C{4,i}(1+par.h(h)+l.cpi:end,h) = log(C{2,i}(startvint+par.h(h)+l.cpi:endvint))-log(C{2,i}(startvint:endvint-par.h(h)-l.cpi));
            else
                C{4,i}(1+par.h(h)+l.ip:end,h) = log(C{2,i}(startvint+par.h(h)+l.ip:endvint))-log(C{2,i}(startvint:endvint-par.h(h)-l.ip));
            end
        end; clear h
        C{4,i} = tseries(startvint:endvint,C{4,i});
        clear startvint endvint
    else
        C{4,i} = C{2,i};
    end
end; clear i
d = struct(C{3:4,:}); % Convert to structure
clear C

% Real-time database
C = cell(4,length(vintages));
C(1:2,:) = [vintages'; struct2cell(DRT)'];
for i=1:length(vintages)
    C{3,i} = lower(C{1,i});
    % Transform levels of indices to common base year
    if  strncmp(C{1,i},'EMP',3)==1 || strncmp(C{1,i},'GDP',3)==1 || strncmp(C{1,i},'UC',2)==1
        C{2,i} = C{2,i};
    elseif strncmp(C{1,i},'DEFL',4)==1 || strncmp(C{1,i},'RENTI',5)==1
        C{2,i} = 100*C{2,i}/mean(C{2,i}(date.baseq));
    else
        C{2,i} = 100*C{2,i}/mean(C{2,i}(date.base));
    end
    % Take h-step log-differences of all time-series except unemployment rate
    if strncmp(C{1,i},'UC',2)==1
        C{4,i} = C{2,i};
    else
        C{4,i} = NaN(size(C{2,i},1),size(par.h,1));
        startvint = get(C{2,i},'start');
        endvint = get(C{2,i},'end');
            if strcmp(C{3,i}(1:3),'cpi') || strcmp(C{3,i}(1:2),'ip')
                for h=1:size(par.h,1)
                    C{4,i}(1+par.h(h)+l.cpi:end,h) = log(C{2,i}(startvint+par.h(h)+l.cpi:endvint))-log(C{2,i}(startvint:endvint-par.h(h)-l.cpi));
                end
            else
                for h=1:size(par.h,1)
                    C{4,i}(1+par.h(h):end,h) = log(C{2,i}(startvint+par.h(h):endvint))-log(C{2,i}(startvint:endvint-par.h(h)));
                end
                C{4,i}(:,1) = [];
            end
        C{4,i} = tseries(startvint:endvint,C{4,i});
    end
end; clear i startvint endvint
DRT = struct(C{1:2,:});
drt = struct(C{3:4,:});
clear C

% Transform quarterly vintages of unemployment to monthly vintages
list = fieldnames(drt);
indc = strfind(list,'uc');
ind = find(not(cellfun('isempty', indc)));
list2 = list(ind,1);
U = cell(2,length(list2));
for i=1:length(list2)
    U{1,i} = list2{i,1};
    U{2,i} = eval(strcat('drt.',list2{i,1}));
end; clear i
UM = cell(2,3*length(list2)-2);
UM2 = cell(2,3*length(list2)-2);
% Set month and year to first vintage (to select vintages accordingly)
countmonth = 11;
countyear = 65;
% Set counter that controls the horizontal movement across VINTAGES!
countvint = 1;
% First observation
date.start = mm(1947,1);
% Set counter that controls the vertical movement across OBSERVATIONS!
countobs = 0.5;
for i=1:3*length(list2)-2
    if countyear<10
        % Add zero in string for two-digit month label
        UM{1,i} = strcat('uc0',num2str(countyear),'m',num2str(countmonth));
    else UM{1,i} = strcat('uc',num2str(countyear),'m',num2str(countmonth));
    end
    % Obtain time-series of vintage at countind
    trash = eval(strcat('drt.',U{1,countvint}));
    % Select observations from that vintage
    UM{2,i} = trash(date.start:date.start+225+i);
    clear trash
    % Add observations from following vintages if observation was de facto
    % available at monthly vintage
    if mod(countobs,1)>0 && i>1
        % Load data from old vintage
        UM{2,i} = UM{2,i-1};
        % Load new observations (max 2!) from following vintage
        trash = eval(strcat('drt.',U{1,countvint+1}));
        UM{2,i}(225+i) = trash(date.start+225+i-1);
        UM{2,i} = [UM{2,i}; NaN];
    end
    clear trash
    % Progress month and year counter
    if countmonth == 12 && countyear==99
        countyear = 0;
        countmonth = 1;
    elseif countmonth>11
        countmonth = 1;
        countyear = countyear+1;
    else countmonth = countmonth+1;
    end
    % Progress observation and vintage counter
    countobs = i/3;
    if mod(countobs,1)==0
        countvint = countvint+1;
    end
    % Convert to time-series
    time = mm(1947,1):mm(1965,11)+i-1;
    UM2{1,i} = UM{1,i};
    UM2{2,i} = tseries(time,UM{2,i});
end
um = struct(UM2{:,:});
for i=1:3*length(list2)-2
    drt.(UM2{1,i}) = um.(UM2{1,i});
end
drt = rmfield(drt,list2);
clear i list list2 countind count countmonth countyear U ind indc UM UM2 um time countobs countvint h

save(strcat('Data\DataOLS',country,'.mat'))