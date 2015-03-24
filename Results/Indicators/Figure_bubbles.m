clear
clc

assetname = 'Stock Prices';
cd(strcat('Indicators\',assetname))

list = dir;
list = list(3:end);
list2 = cell(length(list)-1,1);
for i=1:length(list)-1
    list2{i} = list(i,1).name;
    load(list2{i})
end
clear i list list2

rangeunb = range(p);
startdate = rangeunb(1);
enddate = rangeunb(end);

startyear = dat2str(startdate);
startyear = str2double(startyear{1,1}(1:4));
timeline = startyear+1/12:1/12:startyear+length(p)/12;

assetname = 'Stock Prices';

subplot(3,2,1)
bar(timeline,IndComb3(:))
title(strcat(assetname,': Combi3'),'FontSize',12)
axis('tight')
subplot(3,2,2)
bar(timeline,IndComb4(:))
title(strcat(assetname,': Combi4'),'FontSize',12)
axis('tight')
subplot(3,2,3)
bar(timeline,IndCombPhil2(:))
title(strcat(assetname,': Combi Phillips 2'),'FontSize',12)
axis('tight')
subplot(3,2,4)
bar(timeline,IndCombPhil3(:))
title(strcat(assetname,': Combi Phillips 3'),'FontSize',12)
axis('tight')
subplot(3,2,5)
bar(timeline,[NaN(tau0PSY13-1,1);IndPSY13comb(:)])
title(strcat(assetname,': Combi PSY13'),'FontSize',12)
axis('tight')

cd('..\..\')