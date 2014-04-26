% This function takes the xlsx file 'Longitudinal_RawCounts_ForClass.xlsx'
% and converts it into a 4D numerical matrix.  The 4 dimensions are,
% respectively, UniProtID, mouse ID, colonization state, and location along
% the GI tract.  As part of the process, this function associates each
% dimension's values with a numerical value, e.g. mouse ID==mouse1==1,
% UniProtID==Q8C5B4==1.  These associations are stored in the cell array
% 'axes'.
% The input variable is the filename of the xlsx sheet as a
% string.  The resulting 4 outputs are: 1) The resulting 4D matrix
% 'OverlordMatrix', 2) a Map container'PeptideMap' which returns a matrix
% [mouse ID index, colonization state index, location index] for a specific
% sample of form '{Letter} Peptides'.  The index refers to the numerical
% value associated with a given value as derived from 'axes', e.g. mouse ID
% index = 1for mouse1, as axes{2}{1}=mouse1 (mouse ID is the 2nd dimension and mouse1
% specifically occurs at index 1 of the 2nd dimension), 3) a Map container
% 'LetterMap' which is the inverse of PeptideMap and requires an input of
% the form {mouseID}_{colonization state}_{location}, e.g.
% 'mouse3_RF_ileum', and returns the peptide sample name, e.g. returns 'AE
% Peptides' for the 'mouse3_RF_ileum' input, and 4) the axes cell array as
% described earlier.
function [OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename)
% Rips the 2nd sheet from the xlsx file and generates a numerical array of
% all numerical values, a struct containing all the string values, and a
% struct containing both the numerical and string values.
[num_RawCounts,txt_RawCounts,raw_RawCounts] = xlsread(filename,2);
% Generates the string structs for each colonization state sheet in the
% xlsx file, as peptide samples are associated with a unique colonization
% state.
[num_GF,txt_GF,raw_GF] = xlsread(filename,3,'C1:Q2');
[num_BT,txt_BT,raw_BT] = xlsread(filename,4,'C1:Q2');
[num_RF,txt_RF,raw_RF] = xlsread(filename,5,'C1:Q2');
clear num_GF num_RT num_RF raw_GF raw_BT raw_RF
% For OverlordMatrix, we want the dimensions to be allocated as such:
% [UniProtID,mouseID,colonization state,location].  The signal in each cell
% represents the raw count at the given coordinates,
% e.g. 9001 = [1,1,1,1] means we have a count of 9001 for UniProtID Q8C5B4,
% mouse 3, BT colonized, prox colon]
% To keep things from getting too confusing, I have written in arrays which
% link the set of a given variable to its respective numerical value, e.g.
% axis4key = {cecum,ileum,jejunum,prox colon,stomach} and axis3key =
% {GF,BT,RF}.  Thus, since axis4key{1} = cecum, we know that
% OverlordMatrix(:,:,:,1) contains ALL the data attributed to the cecum.
% Likewise, since axis3key{3} = BT, OverlordMatrix(:,:,3,:) contains ALL
% the data attributed to the BT colonized state.  Likewise,
% OverlordMatrix(:,:,3,1) contains the data attributed to BOTH the cecum
% and the BT colonized state.  This seems a bit cumbersome, but as noted in
% class, it is generally easier to deal with numerical matrices than
% arrays/structs.  Also generates the individual dimensions of the output
% 'axes' struct which is used to convert non-numerical values to associated
% ints.
for ii = 2:1:size(txt_RawCounts,1)
    current_axis1 = txt_RawCounts(ii,2);
    if isnan(current_axis1{1}) == 0 
       axis1key{ii-1} = current_axis1{1};
    else
       current_axis1 = txt_RawCounts(ii,1);
       axis1key{ii-1} = current_axis1{1};
    end
end
axis2key = {'mouse1','mouse2','mouse3'};
axis3key = {'GF','BT','RF'};
% Don't change Stomach capitlization, will break code!
axis4key = {'cecum','ileum','jejunum','prox colon','Stomach'};
% We will now briefly use the Map container to decode the appropriate
% coordinates for every peptide type,  We realize: 1) Each peptide type (A,
% B, C, etc.) is unique, 2) Each peptide type is associated with mouse ID,
% colonization state, and location, 3) we can assign the vector [mouse ID,
% colonization state, location] to each peptide type.  Hence, each peptide
% type functions as a 'key' which returns a 'value', i.e. if Q Peptide,
% then [1,2,1] == [mouse 1, BT, cecum] via our axiskey vectors.  
%%
axes = {axis1key,axis2key,axis3key,axis4key};
PeptideMap = containers.Map();
LetterMap = containers.Map();
for colonizationArray = {txt_BT,txt_GF,txt_RF}
    % First sees which colonization state the current array is associated
    % with and records the colonization state.
    if strcmp(colonizationArray{1}{1,1}, txt_BT{1,1}) == 1;
        % http://stackoverflow.com/questions/8061344/how-to-search-for-a-string-in-cell-array-in-matlab
        % Above is a reference for locating indeces in cell arrays
        colonizationValue = find(ismember(axis3key,'BT'));
        letterColonizationValue = 'BT';
    elseif strcmp(colonizationArray{1}{1,1}, txt_GF{1,1}) == 1;
        colonizationValue = find(ismember(axis3key,'GF'));
        letterColonizationValue = 'GF';
    else
        colonizationValue = find(ismember(axis3key,'RF'));
        letterColonizationValue = 'RF';
    end
    % Associates each peptides sample with a mouse ID and a location and
    % stores data in the PeptideMap.
    for iii = 1:1:size(colonizationArray{1},2)
        % Sets the 'key' for the PeptideMap equal to the sample name, e.g.
        % 'Peptides A'
        key = colonizationArray{1}{1,iii};
        % Examines the value corresponding to the key which contains both
        % the location and the mouse ID
        rawLocationMouse = colonizationArray{1}{2,iii};
        count = length(rawLocationMouse);
        mouseIDValue = '';
        % Since the value containing both the location and mouse ID is
        % writted in the form {location}{mouse ID NUMBER}, e.g. cecum3, we
        % rip off the number at the end and assign it to the mouse ID and
        % assign the remainder to the location, e.g. location = cecum,
        % mouse ID = 3 --> mouse3
        while isstrprop(rawLocationMouse(count),'digit') == 1
            mouseIDValue = strcat(rawLocationMouse(count),mouseIDValue);
            count = count - 1;
            rawLocationMouse = rawLocationMouse(1:count);
        end
        mouseIDValue = strcat('mouse',mouseIDValue);
        locationValue = rawLocationMouse;
        % Takes the non-numerical values for location, mouse ID, and
        % colonization state and assigns each dimension a numerical value
        % by finding the corresponding numerical value in the axes cell array.
        valueVector = [find(ismember(axis2key,mouseIDValue)),colonizationValue,find(ismember(axis4key,locationValue))];
        % Assigns the numerical values to the PeptideMap and records the
        % inverse to the LetterMap
        PeptideMap(key) = valueVector;
        LetterMapKey = strcat(mouseIDValue,'_',letterColonizationValue,'_',locationValue);
        LetterMap(LetterMapKey) = key;
    end
end
%% Generation of the OverlordMatrix
% Iterates throughout the combined numerical/string cell array of the raw
% xlsx data and stores the count value to its respective dimensional
% values.
OverlordMatrix = [];
for ii = 2:1:size(raw_RawCounts,1)
    for iii = 4:1:size(raw_RawCounts,2)
        % Finds the current UniProt ID and converts it to the corresponding
        % numerical value using 'axes' cell array.  Also checks if is is
        % the current protein is a false positive protein as indicated by a
        % "UniProt ID" of NaN.  Uses the full protein name (column 1 in
        % xlsx) instead in this case.
        OverlordMatrix_currentProtein = raw_RawCounts{ii,2};
        if isnan(OverlordMatrix_currentProtein) == 1
           OverlordMatrix_currentProtein = raw_RawCounts{ii,1}; 
        end
        % Returns the current peptide sample name of the current cell and
        % returns the numerical dimensional values associated with the
        % peptide.
        OverlordMatrix_currentPeptide = raw_RawCounts{1, iii};
        OverlordMatrix_currentPeptide = PeptideMap(OverlordMatrix_currentPeptide);
        % Returns the current protein of the current cell and returns the
        % numerical dimensional value associated with the protein.
        OverlordMatrix_currentProtein = find(ismember(axis1key,OverlordMatrix_currentProtein));
        % Assigns the count value of the current cell to the associated
        % numerical dimensional values.
        OverlordMatrix(OverlordMatrix_currentProtein,OverlordMatrix_currentPeptide(1),OverlordMatrix_currentPeptide(2),OverlordMatrix_currentPeptide(3))=raw_RawCounts{ii,iii};
    end
end
end
