% At its highest level, this script first finds which GO IDs occur--since
% there is a theoretical max of 9,999,999 unique GO ID's, it is efficient
% to only find those which occur.  Next, this script assigns a unique numerical
% value to each GO ID, e.g. 1 == GO ID 15.  Lastly, this script finds
% tallies the normalized amount of GO IDs for each GO ID for each sample
% and saves this data to file.  The final output is a 1 x 3013 x 3 x 3 x 5
% matrix, where 3013 corresponds to the number of GO IDs, and 3 3 5
% correspond to the number of samples.  Note that this essentially converts
% from the 876 non-decoy proteins to 3013 GO IDs.  This script uses the
% normOverlord_shannon file.
clear all
close all
clc
load('GOArray')
load('normOverlord_shannon');
load('axes')
GONum = {GOArray{:,3}};
GOOccur = zeros(1,9999999);
min = 9999999;
% Populates a vector to find the occurences of GO ID's; it is known that
% the minimum GO ID value is 15 and the max is like 2e6.  There are a ton
% of empty spaces in between (e.g. GO 27 may not occur), so rather than
% deal with a vector on the order of 1e6 elements, find the values of only
% the occuring GO IDs and sort them into a smaller vector with indeces
% corresponding to GO ID (e.g. index 1 == GO ID 15).  Same reasoning as
% OverlordMatrix and axes.
for ii = 1:1:length(GONum)
    for iii = 1:1:length(GONum{ii})
        GOOccur(GONum{ii}(iii)) = GONum{ii}(iii);
        if GONum{ii}(iii) < min
            min = GONum{ii}(iii);
        end
    end
end
sortedGOOccur = sort(GOOccur);
minindex = find(sortedGOOccur == min);
sortedGOOccur = sortedGOOccur(minindex:length(sortedGOOccur));
clear GOOccur
GOtoIndexConverter = containers.Map();
GOtoIndexConverterStr = containers.Map();
IndextoGOConverter = containers.Map();
IndextoGOConverterStr = containers.Map();
%% Generates Maps to convert between numerical GO ID and index and vice versa
for ii = 1:1:length(sortedGOOccur)
    tempGO = num2str(sortedGOOccur(ii));
    tempIndex = num2str(ii);
    GOtoIndexConverter(tempGO) = ii;
    IndextoGOConverter(tempIndex) = sortedGOOccur(ii);
    zerocount = 7 - length(tempGO);
    tempGO2 = tempGO;
    for j = 1:1:zerocount
        tempGO2 = strcat('0',tempGO2);
    end
    tempGO2 = strcat('GO:',tempGO2);
    IndextoGOConverterStr(tempIndex) = tempGO2;
    GOtoIndexConverterStr(tempGO2) = ii;
end
% We want to find the distribution of GO codes in each sample to see what
% the statics and dynamics are of the GO codes
%%
% Generates a "standard" GO matrix; it is known from the above that there
% are a total of 3013 unique GO codes in this entire dataset.  This chunk
% tallies up the occurence of each GO code for each and every protein in
% the dataset: dimension 1 corresponds to protein, dimension 2 corresponds
% to GO codes.  For example, say protein 1 has GO codes for indeces 3, 500,
% and 1030.  Thus, there will be a value of '1' at proteinGOMatrix(1,3),
% (1,500), and (1,1030) and a value of '0' for all other elements with row
% = 1.
proteinGOMatrix = zeros(size(axes{1},2),length(sortedGOOccur));
for ii = 1:1:size(axes{1},2)
    for iii = 1:1:length(GONum{ii})
        tempString  = num2str(GONum{ii}(iii));
        tempIndex = GOtoIndexConverter(tempString);
        proteinGOMatrix(ii,tempIndex) = proteinGOMatrix(ii,tempIndex) + 1;
    end
end
% Multiplies each GO vector (tallies of GO codes for each protein) by the
% normalized values from the normalized OverlordMatrix.
repproteinGOMat_shannon = repmat(proteinGOMatrix,1,1,3,3,5);

for ii= 1:1:size(normOverlord_shannon,1)
    for iii = 1:1:size(normOverlord_shannon,2)
        for iv = 1:1:size(normOverlord_shannon,3)
            for ivi = 1:1:size(normOverlord_shannon,4)
                repproteinGOMat_shannon(ii,:,iii,iv,ivi) = normOverlord_shannon(ii,iii,iv,ivi) * repproteinGOMat_shannon(ii,:,iii,iv,ivi);
            end
        end
    end
end

% Sums the normalized GO tallies across all proteins, so there is now a
% vector of GO tallies for each sample (sample defined as mouse ID,
% colonization state, and location).
GOenrichMatTemp_shannon = sum(repproteinGOMat_shannon,1);
GOenrichMat_shannon = zeros(size(GOenrichMatTemp_shannon,2),size(GOenrichMatTemp_shannon,3),size(GOenrichMatTemp_shannon,4),size(GOenrichMatTemp_shannon,5));
% Reshapes to be a 4D matrix rather than a 5D

for ii = 1:1:size(GOenrichMatTemp_shannon,3)
    for iii = 1:1:size(GOenrichMatTemp_shannon,4)
        for iv = 1:1:size(GOenrichMatTemp_shannon,5)
            GOenrichMat_shannon(:,ii,iii,iv) = GOenrichMatTemp_shannon(1,:,ii,iii,iv)';
        end
    end
end
save('GOenrichMat_shannon','GOenrichMat_shannon')