%% Hierarchical Clustering
% The following code creates a 2D matrix of proteins vs samples and creates
% a clustergram to visualize hierarchical clustering of the data.

clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = OverlordNormalizer(2); %normalize by row
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
% The following code is still a work in progress.

% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 
clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0005);

% clustergram(all_samples, 'RowPDist', 'spearman');
%%
% The following code is very similar to the clustergram above, but it takes
% a manual approach. It was adapted from Tiffany's code from class.
% Each cluster is graphed on a separate plot to visualize the number of 
% samples in each cluster. 
corrDist = pdist(all_samples);
clusterTree = linkage(corrDist, 'average');
clusters=cluster(clusterTree, 'maxclust', 15); %15 comes from 3 colonization states * 5 locations
 for c = 1:15
     subplot(3,5,c);
     plot(all_samples((clusters == c))');
     axis tight
 end
 suptitle('Hierarchical Clustering of Profiles');

