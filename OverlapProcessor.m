% Retrieves information about overlapping proteins.  Input xlsx must have
% each overlapping sample as a separate column, with each column containing
% the overlapping proteins.  Generates two cell arrays: OverlapGOcodes
% lists the GO codes occurring in each overlap case, while
% OverlapAnnotAndCount lists the corresponding GO code title and counts.
clear all
close all
clc
load axes
load GOArray
load allGODic
[num,txt,raw] = xlsread('RearrangedOverlappingProteins.xlsx');
for i = 1:1:size(txt,2)
    UniProts = {};
    counter = 1;
    while counter <= size(txt,1) && length(txt{counter,i}) > 0
        UniProts{counter} = txt{counter,i};
        counter = counter + 1;
    end
    GOinformation{i} = getGOcodes(UniProts,GOArray,axes,allGODic);
    OverlapGOcodes{i} = keys(GOinformation{i});
    OverlapAnnotAndCount{i} = values(GOinformation{i});
end
save('OverlapGOcodes.mat','OverlapGOcodes')
save('OverlapAnnotAndCount.mat','OverlapAnnotAndCount')