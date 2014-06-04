% Retrieves information about overlapping proteins.  Input xlsx must have
% each overlapping sample as a separate column, with each column containing
% the overlapping proteins.  Generates two cell arrays: OverlapGOcodes
% lists the GO codes occurring in each overlap case, while
% OverlapAnnotAndCount lists the corresponding GO code title and counts.
% OverlapGOcountsvec stores the incidence of each GO code in a more
% accessible matrix format.
clear all
close all
clc
load axes
load GOArray
load allGODic
[num,txt,raw] = xlsread('RearrangedAllProteins.xlsx');
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
OverlapGOcountsvec = [];
for i = 1:1:size(OverlapAnnotAndCount,2)
    for j = 1:1:size(OverlapAnnotAndCount{i},2)
        OverlapGOcountsvec(i,j) = OverlapAnnotAndCount{i}{j}{2};
    end
end
[sortedCounts,sortedI] = sort(OverlapGOcountsvec,2,'descend');
%% Consider the first 5 most prevalent GO codes in case 8 (the strongest overlap)
counter = 1;
for i = 1:1:size(OverlapAnnotAndCount,2)
    if size(OverlapAnnotAndCount{i},2) < 10
        a = 1:1:size(OverlapAnnotAndCount{i},2);
    else
        a = 1:1:10;
    end
    caseAnnot{i,1} = i;
    caseCounts{i,1} = i;
    for j = a
        caseAnnot{i,(j-1)*2+2} = OverlapAnnotAndCount{i}{sortedI(i,j)}{1};
        caseAnnot{i,(j-1)*2+3} = OverlapAnnotAndCount{i}{sortedI(i,j)}{2};
    end
end

fileID = fopen('kmeansGOAnalysisAll10.dat','w');
formatSpec0 = '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t\n\n\n';
header = {'title1-#','title2-#','title3-#','title4-#','title5-#','title6-#','title7-#','title8-#','title9-#','title10-#'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%d %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t %s %d \t\n\n\n';
[nrows, ncols] = size(caseAnnot);
for row = 1:nrows
    fprintf(fileID,formatSpec,caseAnnot{row,:});
end
fclose(fileID);
%% Record the GO codes of the unique proteins
allUniProtMap = containers.Map();
for i = 1:1:size(txt,2)
    UniProts = {};
    counter = 1;
    while counter <= size(txt,1) && length(txt{counter,i}) > 0
        UniProts{counter} = txt{counter,i};
        if isKey(allUniProtMap,UniProts{counter})
            1;
        else
            allUniProtMap(UniProts{counter}) = 1;
        end
        counter = counter + 1;
    end
end
allUniqueGOInfo = containers.Map();
allUniProts = keys(allUniProtMap);
for i = 1:1:length(allUniProts)
    currentGO = getGOcodes({allUniProts{i}},GOArray,axes,allGODic);
    for j = 1:1:length(keys(currentGO))
        currKeys = values(currentGO);
        if isKey(allUniqueGOInfo, currKeys{j}{1})
            value = allUniqueGOInfo(currKeys{j}{1}) + 1;
            allUniqueGOInfo(currKeys{j}{1}) = value;
        else
            allUniqueGOInfo(currKeys{j}{1}) = 1;
        end
    end
end
allUniqueGOs = keys(allUniqueGOInfo);
allUniqueGOval = values(allUniqueGOInfo);
for i = 1:1:length(allUniqueGOval)
    allUniqueGOVal(i) = allUniqueGOval{i};
end
[sortedUniqueVal, valind] = sort(allUniqueGOVal,'descend');
counter = 1;
for i = valind
    allUniqueGO2{counter} = allUniqueGOs{i};
    counter = counter + 1;
end
% allUniqueGOs = allUniqueGOs{valind};
for i = 1:1:length(keys(allUniqueGOInfo))
    allUniqueGOCell{(i-1)*2 + 1} = allUniqueGO2{i};
    allUniqueGOCell{(i-1)*2 + 2} = sortedUniqueVal(i);
end
fileID = fopen('allUniqueGOcounts.dat','w');
formatSpec0 = '%s \t %s\n\n\n';
header = {'title','count'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%s \t %d\n\n\n';
[nrows, ncols] = size(allUniqueGOCell);
for row = 1:nrows
    fprintf(fileID,formatSpec,allUniqueGOCell{row,:});
end
fclose(fileID);

save('GOinformation.mat','GOinformation')
save('OverlapGOcodes.mat','OverlapGOcodes')
save('OverlapAnnotAndCount.mat','OverlapAnnotAndCount')
save('OverlapGOcountsvec.mat','OverlapGOcountsvec')