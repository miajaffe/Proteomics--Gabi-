% Initialization Script
% Run me!
clear all
close all
clc
filename = 'Longitudinal_RawCounts_ForClass_NoDecoy.xlsx';
[OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename);
save('OverlordMatrix.mat')
save('axes');
save('ProteinMap.mat');
save('LetterMap.mat');
% %% Download associated GO-codes for each and every protein
% % Could include functionality to scrape the GO code matrix from multiple
% % computers by dividing labor via partitioning the for loop.
% go_code_array = {};
% for ii = 1:1:size(axes{1},2)
%     fprintf('%d',ii)
%     url = 'https://www.ebi.ac.uk/interpro/protein/';
%     url = strcat(url,axes{1}{ii});
%     urldata = urlread(url);
%     proteinindeces1 = findstr(urldata,'<p>');
%     proteinindeces2 = findstr(urldata,'</p>');
%     % The ones after the 3rd index are irrelevant
%     if length(proteinindeces1) >= 3
%         proteinindeces1 = proteinindeces1(1:3);
%         proteinindeces2 = proteinindeces2(1:3);
%         proteinindeces = findstr(urldata,'GO:');
%         proteinindeces(length(proteinindeces)+1) = proteinindeces(length(proteinindeces)) + 50;
%         biocount = 0;
%         molcount = 0;
%         cellcount = 0;
%         % Assume GO Code ID is 7 digits long, after the 'GO:' and that the
%         % first 3 entries are for a sample GO code which do not chnage over
%         % time and that each GO code is mentioned twice.
%         for j = 5:2:length(proteinindeces)
%             currStr = urldata(proteinindeces(j):proteinindeces(j+1));
%             descriploc1 = findstr(currStr,'>');
%             descriploc2 = findstr(currStr,'<');
%             descriploc1 = descriploc1(1);
%             if length(descriploc2) > 1
%                 descriploc2 = descriploc2(2);
%             else
%                 descriploc2 = length(currStr);
%             end
%             if proteinindeces(j) < proteinindeces2(1)
%                 biocount = biocount + 1;
%                 go_bio{ii}{biocount}{1} = urldata(proteinindeces(j)+3:proteinindeces(j)+9);
%                 go_bio{ii}{biocount}{2} = currStr(descriploc1+2:descriploc2-1);
%             elseif proteinindeces(j) > proteinindeces2(1) && proteinindeces(j) < proteinindeces2(2)
%                 molcount = molcount + 1;
%                 go_mol{ii}{molcount}{1} = urldata(proteinindeces(j)+3:proteinindeces(j)+9);
%                 go_mol{ii}{molcount}{2} = currStr(descriploc1+2:descriploc2-1);
%             elseif proteinindeces(j) > proteinindeces2(2) && proteinindeces(j) < proteinindeces2(3)
%                 cellcount = cellcount + 1;
%                 go_cell{ii}{cellcount}{1} = urldata(proteinindeces(j)+3:proteinindeces(j)+9);
%                 go_cell{ii}{cellcount}{2} = currStr(descriploc1+2:descriploc2-1);
%             end
%         end
%     end
%     %     count = 0;
%     %     for iii = 4:2:length(proteinindeces)
%     %         count = count + 1;
%     %         go_code_array{ii}{count} = urldata(proteinindeces(iii)+3:proteinindeces(iii)+9);
%     %     end
% end
% save('go_bio.mat','go_mol.mat','go_cell.mat');
% % As a note, as of 4/29/14, the final protein, Q8CAQ8 returns no associated
% % GO ID's, so the go_code_array, while only having 891 cells filled
% % relative to the 892 proteins, is truncated and not rearranged, which
% % makes processing easier.
% %% Convert GO code ID strings to numbers
% for ii = 1:1:length(go_code_array)
%     for iii = 1:1:length(go_code_array{ii})
%         temp_vector(iii) = str2num(go_code_array{ii}{iii});
%     end
%     go_code_array_num{ii} = temp_vector;
% end
% save('go_codes_num.mat','go_code_array_num');
