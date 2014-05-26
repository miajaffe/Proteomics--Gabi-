%% Hierarchical Clustering of GO codes
clear all
close all hidden
clc
load('axes140523.mat');
load('normOverlordFinal_140523.mat');
load('GOenrichMat'); % 3013 * 3 * 3 * 5 4d matrix NOTE: This isn't the most updated version!
all_samples = [];
all_labels = {};
normOverlord = normOverlordFinal; 

% This loops through all mouse #s, colonization states, and locations to
% create a 2D matrix where each row is a gocode_id and each column is a
% different sample (sample A, B, etc.). The matrix, all_samples, contains the abundance
% of each GO code at each sample.  The all_labels matrix is a 1x45 matrix
% that contains labels for each sample in the order that they are in teh
% all_samples matrix. 
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            all_samples = [all_samples GOenrichMat(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end
    end
end

%simple stats - needed to determine cutoff of clustergram
num_samples = length(GOenrichMat(:,1,1,1)) * 3 * 3 * 5;
reshaped_samples = reshape(all_samples, num_samples,1);
percentages = prctile(reshaped_samples, [25 50 75 90 95]);
% 50% = 0 75% = 0.0006  90% = 0.0036   95% = 0.0084
percentages2 = prctile(reshaped_samples, [81 82 83 84 85]);
% 81% = 0.0011 83% = 0.0014 84% = 0.0016

%This clustergram is capturing between 83-84% of the data. In other words
%16-17% of the data is above the 0.0015 cutoff.
cluster = clustergram(all_samples, 'ColumnLabels', all_labels', 'DisplayRange', 0.0015, 'Symmetric', 'true', 'Colormap', winter);
