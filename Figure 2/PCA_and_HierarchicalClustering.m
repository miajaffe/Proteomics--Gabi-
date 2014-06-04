
% Figure 2

% Load Initialized Data
clear all;
close all hidden;
clc;
load('axes.mat');
% load('LetterMap.mat');
load('../Initialization/normOverlordFinal.mat');
%load('ProteinMap.mat');
% load('MusProt.mat');
all_samples = [];
all_labels = {};

% Loops through all mouse replicates, colonization states, and locations to
% create a 2D matrix where each row is a protein_id and each column is a
% different sample (sample A, B, etc.). The matrix, all_samples, contains the abundance
% of each protein at each sample.  The all_labels matrix is a 1x45 matrix
% that contains labels for each sample in the order that they are in teh
% all_samples matrix. The all_labels matrix and the proteins matrix can be 
% used for labeling the axes of any plots.
for mouse_num = 1:3;
    for colonization = 1:3;
        for loc = 1:5;
            all_samples = [all_samples normOverlordFinal(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end;
    end;
end;
            

%% PCA 
% By Daniel Sprockett
% step 1: Calculate the average over all vectors
avg = mean(all_samples);

%%
number_of_proteins = size(all_samples,1);
number_of_samples = size(all_samples,2);

% step 2: Subtract the average from all vectors
diff_avg = all_samples-repmat(avg, number_of_proteins, 1);

%%
% step 3: Calculate the covariance_matrixariance matrix

covariance_matrix = zeros(number_of_samples, number_of_samples);

for i = 1:number_of_samples;
    for j = 1:number_of_samples;
        covariance_matrix(i,j) = 1/(number_of_proteins - 1) * ...
            sum(diff_avg(:,i).*diff_avg(:,j));
    end;
end;
figure;
imagesc(covariance_matrix);
colorbar;

%%
% step 4: calculate eigs of covariance_matrix

[eigenvectors,eigenvalues,explained] = pcacov(covariance_matrix);


%%

% Create Color Maps

for i = 1: numel(all_labels);
    split_strings(i,:) = strsplit(all_labels{1,i}, '_');
end

% Gut Location

color_by_GutRegion = zeros(45, 3);

 for i = 1: numel(all_labels);

    if strcmp(split_strings(i,3), 'cecum'); % Color cecum Red
    color_by_GutRegion(i,1) = 1; color_by_GutRegion(i,2) = 0; color_by_GutRegion(i,3) = 0;

    elseif strcmp(split_strings(i,3), 'ileum'); % Color ileum Green
    color_by_GutRegion(i,1) = 0; color_by_GutRegion(i,2) = 1; color_by_GutRegion(i,3) = 0;

    elseif strcmp(split_strings(i,3), 'jejunum'); % Color jejunum Blue
    color_by_GutRegion(i,1) = 0; color_by_GutRegion(i,2) = 0; color_by_GutRegion(i,3) = 1;

    elseif strcmp(split_strings(i,3), 'prox colon'); % Color prox colon cyan
    color_by_GutRegion(i,1) = 0; color_by_GutRegion(i,2) = 1; color_by_GutRegion(i,3) = 1;

    else strcmp(split_strings(i,3), 'Stomach'); % Color Stomach magenta
    color_by_GutRegion(i,1) = 1; color_by_GutRegion(i,2) = 0; color_by_GutRegion(i,3) = 1;
    
    end;
 end;
 
% Colonization State

color_by_ColonizationState = zeros(45, 3);

 for i = 1: numel(all_labels);

    if strcmp(split_strings(i,2), 'BT'); % Color BT mice Green
    color_by_ColonizationState(i,1) = .2; color_by_ColonizationState(i,2) = 1; color_by_ColonizationState(i,3) = .2;

    elseif strcmp(split_strings(i,2), 'RF'); % Color RF mice Red
    color_by_ColonizationState(i,1) = 1; color_by_ColonizationState(i,2) = 0; color_by_ColonizationState(i,3) = .3;

    else strcmp(split_strings(i,2), 'GF'); % Color Germ Free mice Blue 
    color_by_ColonizationState(i,1) = .3; color_by_ColonizationState(i,2) = .5; color_by_ColonizationState(i,3) = 1;
    
    end
 end


%%

% Plot PC1 vs PC2 vs PC3 for Gut Region

figure
scatter3(eigenvectors(:,1), eigenvectors(:,2), eigenvectors(:,3), 100, color_by_GutRegion, 'filled');
hold on;
title('Principal Component Analysis: Gut Location');
xlabel(sprintf('PC1 = %.0f%%', explained(1,1)));
ylabel(sprintf('PC2 = %.0f%%', explained(2,1)));
zlabel(sprintf('PC3 = %.0f%%', explained(3,1)));
[~, I] = unique(split_strings(:,3)); 
I = length(split_strings(:,3)) - I(length(I):-1:1); 
p = findobj(gca,'Type','Patch');
legend(p(I), 'stomach', 'prox colon', 'jejunum', 'ileum', 'cecum');
hold off;


%%

% Plot PC1 vs PC2 for Gut Region

figure
scatter(eigenvectors(:,1), eigenvectors(:,2), 100, color_by_GutRegion, 'filled');
hold on;
%title('Principal Component Analysis: Gut Location');
xlabel(sprintf('PC1 = %.0f%%', explained(1,1)));
ylabel(sprintf('PC2 = %.0f%%', explained(2,1)));
[~, I] = unique(split_strings(:,3)); 
I = length(split_strings(:,3)) - I(length(I):-1:1); 
p = findobj(gca,'Type','Patch');
legend(p(I), 'stomach', 'prox colon', 'jejunum', 'ileum', 'cecum', 'Location','SouthWest');
hold off;

%%
% Plot PC2 vs PC3 for Gut Region

figure
scatter(eigenvectors(:,2), eigenvectors(:,3), 100, color_by_GutRegion, 'filled');
hold on;
%title('Principal Component Analysis: Gut Location');
xlabel(sprintf('PC2 = %.0f%%', explained(2,1)));
ylabel(sprintf('PC3 = %.0f%%', explained(3,1)));
[~, I] = unique(split_strings(:,3)); 
I = length(split_strings(:,3)) - I(length(I):-1:1); 
p = findobj(gca,'Type','Patch');
%legend(p(I), 'stomach', 'prox colon', 'jejunum', 'ileum', 'cecum');
hold off;

%%

% Plot PC1 vs PC2 vs PC3 for Colonization State

figure;
scatter3(eigenvectors(:,1), eigenvectors(:,2), eigenvectors(:,3), 100, color_by_ColonizationState, 'filled');
hold on;
title('Principal Component Analysis: Colonization State');
xlabel(sprintf('PC1 = %.0f%%', explained(1,1)));
ylabel(sprintf('PC2 = %.0f%%', explained(2,1)));
zlabel(sprintf('PC3 = %.0f%%', explained(3,1)));
[~, I] = unique(split_strings(:,2)); 
I = length(split_strings(:,2)) - I(length(I):-1:1); 
p = findobj(gca,'Type','Patch');
legend(p(I), 'BT', 'RF', 'GF');
hold off;


%% 
% According to the princomp help files, the loadings are the COEFF file.
% However, but I don't understand how to interpret them. Also, princomp is
% old and they don't recommend using it. 

[COEFF,SCORE,latent,tsquare] = princomp(covariance_matrix);


%%
% PCA Visulization Tool
% Note: I don't think this is very useful for us, but MATLAB has it
 mapcaplot(all_samples);
 

 
  %% Hierarchical Clustering
% The following code creates a 2D matrix of proteins vs samples and creates
% a clustergram to visualize hierarchical clustering of the data.

 
%%
% By Emily Alsentzer
% The following code is still a work in progress.

% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 

load('axes.mat');
proteins = axes{1}'; %generates list of protein ids

clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0015, 'Symmetric', 'true', 'colormap', 'winter');
% colorbar
% hold on
% addTitle('Hierarchical Clustering of Proteins');
% addXLabel('Samples');
% addYLabel('Proteins');
% hold off
 
% note: add 'RowLabels', proteins,  when the axes file is fixed
% clustergram(all_samples, 'RowPDist', 'spearman');

%%
% The following code is very similar to the clustergram above, but it takes
% a manual approach. It was adapted from Tiffany's code from class.
% Each cluster is graphed on a separate plot to visualize the number of 
% samples in each cluster. 
corrDist = pdist(all_samples);
clusterTree = linkage(corrDist, 'average');
clusters=cluster(clusterTree, 'maxclust', 15); %15 comes from 3 colonization states * 5 locations
 for c = 1:15;
     subplot(3,5,c);
     plot(all_samples((clusters == c))');
     axis tight;
 end
 suptitle('Hierarchical Clustering of Profiles');



