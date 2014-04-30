clear all
close all hidden
clc
% Run this script only if the .mat's derived from Initialization_Script.m
% are in the current folder!
load('axes.mat');
load('go_codes_final.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('ProteinMap.mat');
load('go_codes_num.mat');
% Assume that the first GO ID in every cell corresponds to a particular
% biological process (this should be worked into later versions of the
% Initialization_Script.m so that the biological, molecular, and cellular 
% components are kept separate in the go_code_array
%% Generate a vector of the first biological GO ID (crude, I know)
first_biological_GO = zeros(length(go_code_array_num),1);
for ii = 1:1:length(go_code_array_num)
    first_biological_GO(ii) = go_code_array_num{ii}(1);
end
all_GO_codes = [go_code_array_num{:}];
%% What are the most common biological functions of proteins?
close all
% Condense all proteins which share their biological GO ID into separate
% bins for each GO ID to get an estimate of which functions are commonest
binranges = [1:1:max(first_biological_GO)];
[bincounts,ind] = histc(first_biological_GO,binranges);
[sorted_bins,sorted_indeces] = sort(bincounts,'descend');
non_zero_indeces = find(sorted_bins);
% Truncate sorted_bins and sorted_indeces to remove zero-occuring GO ID's
sorted_bins = sorted_bins(non_zero_indeces);
sorted_indeces = sorted_indeces(non_zero_indeces);
% cited: http://stackoverflow.com/questions/6397234/matlab-how-to-use-strings-instead-of-numbers-in-bar-figure
for ii = 1:1:length(sorted_indeces)
    bar_string{ii} = num2str(sorted_indeces(ii));
end
figure
bar(sorted_bins)
set(gca, 'XTickLabel',bar_string,'XTick',1:length(bar_string))
ylabel('Count of occurences')
xlabel('GO IDs')
title('The Unique Biol. GO IDs in the Proteins identified')
% Yikes!  Let's constrain ourselves to the 5 most prevalent GO ID's for
% starters #3many5me #gg_no_re
figure
five_most_bins = sorted_bins(1:5);
five_most_indeces = sorted_indeces(1:5);
for ii = 1:1:5
    bar_string_five{ii} = num2str(sorted_indeces(ii));
end
bar(five_most_bins)
set(gca, 'XTickLabel',bar_string_five,'XTick',1:length(bar_string_five))
ylabel('Count of occurences')
xlabel('GO IDs')
title('The 5 Most Unique Biol. GO IDs in the Proteins identified')
% This is all fine and dandy, but there may not necessarily be a
% correlation between the proteins listed and the peptides identified.  Idk
% yet, man, idk.  ONLY ONE WAY TO FIND OUT
% BEWARE!  I CANNOT STRESS ENOUGH THAT THIS IS USING THE NON-NORMALIZED
% DATA, THIS IS MOSTLY TO GIVE A ROUGH IDEA OF WHAT THE PLAYING FIELD LOOKS
% LIKE AND ESTABLISH A FRAMEWORK FOR THE FUTURE!

%% Reformat this section to function to find commonest biol. function for given conditions
% Count the total number of occurences for each gene, disregarding location
% and sample and all that jazz
cumulative_counts = sum(sum(sum(OverlordMatrix,4),3),2);

% We need to pool together proteins with equivalent GO ID's to avoid
% redundancy downstream.  Sort the GO ID's in descending order
[sorted_biological_GO,sorted_biological_indeces] = sort(first_biological_GO);
% Sort the cumulative counts to correspond to the original GO ID
sorted_cumulative_counts = cumulative_counts(sorted_biological_indeces);
% We now want to make a new matrix which stores only unique GO ID's and the
% corresponding cumulative count for each GO ID
unique_cum_matrix = [sorted_biological_GO(1),sorted_cumulative_counts(1)];
counter = 1;
for ii = 2:1:length(sorted_cumulative_counts)
    if sorted_biological_GO(ii) == sorted_biological_GO(ii-1)
       unique_cum_matrix(counter,2) = unique_cum_matrix(counter,2) + sorted_cumulative_counts(ii);
    else
        counter = counter + 1;
        unique_cum_matrix(counter,1) = sorted_biological_GO(ii);
    end
end
% Sort unique_cum_matrix according to count (column 2)
[sorted_cum_counts,sorted_cum_indeces] = sort(unique_cum_matrix(:,2),'descend');
% Sort corresponding GO ID
sorted_cum_GO_temp = unique_cum_matrix(:,1);
sorted_cum_GO = sorted_cum_GO_temp(sorted_cum_indeces);
% As before, let's look at the top 5 proteins
cum_biological_GO = sorted_cum_GO(1:5);
figure
for ii = 1:1:5
    cum_bar_string_five{ii} = num2str(cum_biological_GO(ii));
end
bar(sorted_cum_counts(1:5))
set(gca, 'XTickLabel',cum_bar_string_five,'XTick',1:length(cum_bar_string_five))
ylabel('Count of occurences')
xlabel('GO IDs')
title('The 5 Most Unique Biol. GO IDs in the Proteins sampled')
%% What the hell, might as well see what happens at different gut locations
figure
% Stomach = 5, jejunum = 3, ileum = 2, cecum = 1, prox colon = 4
GI_vector = [5,3,2,1,4];
for jj = 1:1:5
% Consider only elements pertaining to element 5 of the 4th dimension
GI_OverlordMatrix = OverlordMatrix(:,:,:,GI_vector(jj));
cumulative_counts = sum(sum(sum(GI_OverlordMatrix,4),3),2);

% We need to pool together proteins with equivalent GO ID's to avoid
% redundancy downstream.  Sort the GO ID's in descending order
[sorted_biological_GO,sorted_biological_indeces] = sort(first_biological_GO);
% Sort the cumulative counts to correspond to the original GO ID
sorted_cumulative_counts = cumulative_counts(sorted_biological_indeces);
% We now want to make a new matrix which stores only unique GO ID's and the
% corresponding cumulative count for each GO ID
unique_cum_matrix = [sorted_biological_GO(1),sorted_cumulative_counts(1)];
counter = 1;
for ii = 2:1:length(sorted_cumulative_counts)
    if sorted_biological_GO(ii) == sorted_biological_GO(ii-1)
       unique_cum_matrix(counter,2) = unique_cum_matrix(counter,2) + sorted_cumulative_counts(ii);
    else
        counter = counter + 1;
        unique_cum_matrix(counter,1) = sorted_biological_GO(ii);
    end
end
% Sort unique_cum_matrix according to count (column 2)
[sorted_cum_counts,sorted_cum_indeces] = sort(unique_cum_matrix(:,2),'descend');
% Sort corresponding GO ID
sorted_cum_GO_temp = unique_cum_matrix(:,1);
sorted_cum_GO = sorted_cum_GO_temp(sorted_cum_indeces);
% As before, let's look at the top 5 proteins
cum_biological_GO = sorted_cum_GO(1:5);
% figure
for ii = 1:1:5
    cum_bar_string_five{ii} = num2str(cum_biological_GO(ii));
end
subplot(3,2,jj)
bar(sorted_cum_counts(1:5))
set(gca, 'XTickLabel',cum_bar_string_five,'XTick',1:length(cum_bar_string_five))
ylabel('Count of occurences')
xlabel('GO IDs')
% Stomach = 5, jejunum = 3, ileum = 2, cecum = 1, prox colon = 4
titletemp = {' stomach',' jejunum',' ileum',' cecum',' prox colon'};
titlestr = strcat('5 commonest GO ID at ' , titletemp{jj});
title(titlestr)

end