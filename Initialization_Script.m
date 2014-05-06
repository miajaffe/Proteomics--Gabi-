% Initialization Script
% Run me!
clear all
close all
clc
filename = 'Longitudinal_RawCounts_ForClass.xlsx';
[OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename);
%% Download associated GO-codes for each and every protein
% Could include functionality to scrape the GO code matrix from multiple
% computers by dividing labor via partitioning the for loop.
go_code_array = {};
for ii = 1:1:size(axes{1},2)
    fprintf('%d',ii)
    url = 'https://www.ebi.ac.uk/interpro/protein/';
    url = strcat(url,axes{1}{ii});
    urldata = urlread(url);
    proteinindeces = findstr(urldata,'GO:');
    % Assume GO Code ID is 7 digits long, after the 'GO:' and that the
    % first 3 entries are for a sample GO code which do not chnage over
    % time and that each GO code is mentioned twice.  Pesumably, urlread2
    % isfaster than urlread (not official matlab script).
    count = 0;
    for iii = 4:2:length(proteinindeces)
        count = count + 1;
        go_code_array{ii}{count} = urldata(proteinindeces(iii)+3:proteinindeces(iii)+9);
    end
end
save('go_codes.mat','go_code_array');
% As a note, as of 4/29/14, the final protein, Q8CAQ8 returns no associated
% GO ID's, so the go_code_array, while only having 891 cells filled
% relative to the 892 proteins, is truncated and not rearranged, which
% makes processing easier.
%% Convert GO code ID strings to numbers
for ii = 1:1:length(go_code_array)
    for iii = 1:1:length(go_code_array{ii})
        temp_vector(iii) = str2num(go_code_array{ii}{iii});
    end
    go_code_array_num{ii} = temp_vector;
end
save('go_codes_num.mat','go_code_array_num');
