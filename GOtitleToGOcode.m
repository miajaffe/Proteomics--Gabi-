% Converts GO title to GO code and index
% Keys are titles, e.g. 'zymogengranule', and values are a 1x3 cell
% array--for this key, the values are: {'GO:0042589', 1977, 4-D double},
% where 1977 is the GO index in GOenrich mat and the 4-D double is the
% 1 x 3 x 3 x 5 matrix containing the normalized abundances across all
% samples for the given GO code.
clear all
close all
clc
load GOenrichMat
load GOenrichMat_shannon
load axes
load GOtoIndexConverterStr
load IndextoGOConverterStr
load allGODic
temp = values(allGODic);
temp2 = keys(allGODic);
titleToGO = containers.Map();
for i = 1:1:size(temp,2)
    value1 = temp2{i};
    value2 = GOtoIndexConverterStr(value1);
    value3 = GOenrichMat(value2,:,:,:);
    titleToGO(temp{i}{1}) = {value1,value2,value3};
end
save('titleToGO.mat','titleToGO')