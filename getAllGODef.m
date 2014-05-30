% GET ALL THE GO CODE TITLES.
clear all
clc
load GOtoIndexConverter
allnumGO = keys(GOtoIndexConverter);
% Convert numerical GO values to searchable strings
allstrGO = {};
for i = 1:1:length(allnumGO)
    currStr = num2str(allnumGO{i});
    front = 'GO:';
    zeroslen = 7 - length(currStr);
    for j = 1:1:zeroslen
        front = strcat(front,'0');
    end
    allstrGO{i} = strcat(front,currStr);
end

% Store annotation and definition in a Map with GO ID as key
allGODic = containers.Map();
for i = 1:1:size(allstrGO,2)
    key = allstrGO{i};
    url = strcat('http://www.ebi.ac.uk/QuickGO/GTerm?id=',key);
    tempData = urlread(url);
    startI = strfind(tempData,'<td class="label">Definition') + 40;
    counter = startI;
    definition = '';
    while tempData(counter) ~= '<';
        definition = strcat(definition,tempData(counter));
        counter = counter + 1;
    end
%     endI = strfind(tempData,'.    </td>  </tr>        <tr>  <td class="label">');
%     roughData = tempData(startI:endI(2));
%     polishedData = roughData(49:length(roughData)-5);
%     definition = polishedData;
    startI = strfind(tempData,'<title>') + 18;
    annot = '';
    counter = startI;
    while tempData(counter) ~= '<';
        annot = strcat(annot,tempData(counter));
        counter = counter + 1;
    end
    allGODic(key) = {annot,definition};
    fprintf('%d\n',i)
end
save('allGODic.mat','allGODic')