%% Redo hierarchical with Mia's proteins only
clear all
close all hidden
clc
load('axes140523.mat');
load('normOverlordFinal_140523.mat');
load('mia_clustering_work_with_final_norm.mat');

normOverlord = normOverlordFinal; 
%%
% In the mat file "mia_clustering_work_with_final_norm.mat", the following 
% three variables contain the indices relative in the OverlordMatrixFinal 
% of the proteins contained in the kmeans clusters:
% 'i_bt_hi'
% 'i_cv_hi'
% 'i_gf_hi'
A = [1:843]

bt_indices = i_bt_hi;
gf_indices = i_gf_hi;
conv_indices = i_cv_hi;
total_indices = [bt_indices; gf_indices; conv_indices];
final_indices = unique(total_indices);
%% if the indices are the high abundant proteins
reduced_overlord = normOverlord(final_indices,:,:,:);

%%
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
            all_samples = [all_samples reduced_overlord(:,mouse_num, colonization, loc)];
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

reshape_size = length(all_samples(:,1)) * length(all_samples(1,:));
reshaped_samples = reshape(all_samples, reshape_size,1);
percentages = prctile(reshaped_samples, [25 50 75 90 95])
% % % 50% = 0 75% = 0.0002  90% = 0.0015   95% = 0.0043
percentages2 = prctile(reshaped_samples, [81 82 83 84 85])
% 

           
%%
% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 

% 10% of the data is above the cutoff.
cluster = clustergram(all_samples, 'ColumnLabels', all_labels', 'DisplayRange', 0.0271, 'Symmetric', 'true', 'Colormap', winter);
%cluster = clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0015, 'Symmetric', 'true', 'Colormap', winter);




