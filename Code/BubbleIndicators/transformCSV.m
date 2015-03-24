clear
clc

cd('..\Daten\Non-Revised\')

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

list = list2(5);

for j=1:length(list)
    namein = list{j}(1:end-4);
    nameout = strcat(namein,'_new');

    fid0 = fopen(strcat(namein,'.csv'),'r');
    fid1 = fopen(strcat(nameout,'.csv'),'w');
    fwrite(fid1,strrep(char(fread(fid0))',';',','));
    fclose(fid0);
    fclose(fid1);
end

cd('..\..\Code')