% This function takes a list of UniProt IDs of interest as well as the axes
% and GOArray matrices.  The function sequentially finds the GO codes of
% each UniProt ID and stores the information in a Map structure.  The keys
% are the GO IDs found and the values are a 1 x 2 cell array, containing in
% the first cell the annotation and the total incidence of said GO ID
% within the list of UniProt IDs in the 2nd cell.  For example: key =
% GO:99999999, value = {'unknown function',2}, given that two proteins have
% the GO ID GO:99999999.  Also requires allGODic.

%Must load GOArray, axes, allGODic
function [GOinformation] = getGOcodes(UniProts,GOArray,axes,allGODic)
GOinformation = containers.Map();
for i = 1:1:length(UniProts)
    currentUniProt = UniProts{i};
    currProtIndex = find(ismember(axes{1},currentUniProt));
    currGOIDs = GOArray{currProtIndex,1};
    for j = 1:1:length(currGOIDs)
        if isKey(GOinformation,currGOIDs{j})
            % Value is a 1x2 cell array with tally in the 2nd cell
            value = GOinformation(currGOIDs{j});
            value{2} = value{2} + 1;
            GOinformation(currGOIDs{j}) = {value{1},value{2}};
        else
            annot = allGODic(currGOIDs{j});
            GOinformation(currGOIDs{j}) = {annot{1},1};
        end
    end
end
end