clear
clc

cd('..\Daten\Realtime Daten')

list = dir;
numvars = length(list)-2;
list2 = cell(numvars,1);
list3 = cell(numvars,1);
for i=1:numvars
    list2{i} = list(i+2,1).name;
    list3{i} = list2{i}(1:end-4);
end; clear i

list = list3;
clear list2 list3

for i=267:numvars
    [data,textdata] = readtable(list{i});
    start = textdata(2,1);
    starty = start{1}(1:4);
    startqm = start{1}(end-1:end);
    if strcmp(startqm(1),'0')
        startqm = startqm(2);
    end
    assignin('base',list{i},[textdata(end-size(data,1)+1:end,1),data]);
end