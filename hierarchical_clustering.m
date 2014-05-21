%% Hierarchical Clustering
% The following code creates a 2D matrix of proteins vs samples and creates
% a clustergram to visualize hierarchical clustering of the data.

clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('normOverlord2.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = normOverlord2; 
all_samples = [];
all_labels = {};
proteins = axes{1}'; %generates list of protein ids

% This loops through all mouse #s, colonization states, and locations to
% create a 2D matrix where each row is a protein_id and each column is a
% different sample (sample A, B, etc.). The matrix, all_samples, contains the abundance
% of each protein at each sample.  The all_labels matrix is a 1x45 matrix
% that contains labels for each sample in the order that they are in teh
% all_samples matrix. The all_labels matrix and the proteins matrix can be 
% used for labeling the axes of any plots.
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            all_samples = [all_samples normOverlord(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end
    end
end
%%
%Simple statistics
ranges = range(all_samples);
max_range = max(ranges) %range = 0.2343 
boxplot(all_samples);
title('Distribution of normalized abundance per sample')
xlabel('Sample')
ylabel('Normalized Protein Abundance')
median = median(all_samples) % median is 0 for all except col 32 where median = 0.1402
percentages = prctile(all_samples, [25 50 75 90 95]); %percentages for each col

reshaped_samples = reshape(all_samples, 39420,1);
percentages = prctile(reshaped_samples, [25 50 75 90 95])
% 50% = 0 75% = 0.0002  90% = 0.0015   95% = 0.0043
percentages2 = prctile(reshaped_samples, [81 82 83 84 85])


%%Should we remove proteins with 0 values?!

           
%%
% The following code is still a work in progress.

% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 

% 10% of the data is above the cutoff.
cluster = clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0015, 'Symmetric', 'true', 'Colormap', winter);

% clustergram(all_samples, 'RowPDist', 'spearman');
%%
% The following code is very similar to the clustergram above, but it takes
% a manual approach. It was adapted from Tiffany's code from class.
% Each cluster is graphed on a separate plot to visualize the number of 
% samples in each cluster. 
corrDist = pdist(all_samples);
clusterTree = linkage(corrDist, 'average');
clusters=cluster(clusterTree, 'maxclust', 5); 
 for c = 1:5
     subplot(5,1,c);
     plot(all_samples((clusters == c))')
     axis tight
 end
 suptitle('Hierarchical Clustering of Profiles');

