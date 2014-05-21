
% By Emily Alsentzer
clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('ProteinMap.mat');
load('GOenrichMat.mat');
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
    end
end
figure
imagesc(covariance_matrix);
colorbar;

%%
% step 4: calculate eigs of covariance_matrix

[eigenvectors,eigenvalues,explained] = pcacov(covariance_matrix);

%%
% Create Color Maps

% Gut Location

for i = 1: numel(all_labels);
    split_strings(i,:) = strsplit(all_labels{1,i}, '_');
end

color_by_GutRegion_strings = split_strings(:,3);

% color_by_GutRegion = zeros(45, 3);

cecum_color = strcmp(color_by_GutRegion_strings, 'cecum');

if cecum_color == 1;
   pass = 1;
end

% elseif
    'ileum'
    'jejunum'
    'prox colon'
    'Stomach'

color_by_ColonizationState = split_strings(:,2);

% Plot PC1 vs PC2 vs PC3

scatter3(eigenvectors(:,1), eigenvectors(:,2), eigenvectors(:,3), 100, color_by_GutRegion, 'filled');

%% 
% step 5: generate files needed for a PCA plot

eigenvector_matrix = eigenvectors;
add_labels = {'eigvals'; '% variation explained'};
add_title = {'pc vector number'};
pc_vector_labels = all_labels';
pc_vector_number = vertcat(add_title,pc_vector_labels, add_labels);

pc_labels = (1:size(all_labels,2));
eigvals = eigenvalues';
variation_explained = explained';

eigenvector_matrix_labeled = vertcat(pc_labels, eigenvector_matrix, eigvals, variation_explained);
eigenvector_table = table(pc_vector_number, eigenvector_matrix_labeled);

writetable(eigenvector_table,'QIIME_proteomics_pc.txt','Delimiter','\t');

% Then open in microsoft excel and delete the first line.

%%
% plot data with the eigenvectors 
figure;
scatter3(all_samples(:,1),all_samples(:,2),...
    all_samples(:,3),100,[1 0 0],'filled');
ev1 = eigenvectors(:,1)*eigenvalues(1,1);
ev2 = eigenvectors(:,2)*eigenvalues(2,2);
ev3 = eigenvectors(:,3)*eigenvalues(3,3);
hold on;
%plot3([avg(1)-ev1(1) avg(1)+ev1(1)],[avg(2)-ev1(2) avg(2)+ev1(2)], ...
 %    [avg(3)-ev1(3) avg(3)+ev1(3)],'Color',[0 0 1],'LineWidth',3);
% plot3([avg(1)-ev2(1) avg(1)+ev2(1)],[avg(2)-ev2(2) avg(2)+ev2(2)], ...
%     [avg(3)-ev2(3) avg(3)+ev2(3)],'Color',[0 0 0],'LineWidth',3);
%plot3([avg(1)-ev3(1) avg(1)+ev3(1)],[avg(2)-ev3(2) avg(2)+ev3(2)], ...
 %    [avg(3)-ev3(3) avg(3)+ev3(3)],'Color',[0 1 1],'LineWidth',3);
hold off;
axis equal

%%
 mapcaplot(all_samples)
 
 %%
 
 
  %% Hierarchical Clustering
% The following code creates a 2D matrix of proteins vs samples and creates
% a clustergram to visualize hierarchical clustering of the data.

 
%%
% By Emily Alsentzer
% The following code is still a work in progress.

% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 
clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0005, 'colormap', 'winter');

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


%%


raw_data = xlsread('Longitudinal_RawCounts_ForClass.xlsx', 2);

number_of_proteins = size(raw_data,1);
number_of_samples = size(raw_data,2);

raw_data_avg = mean(raw_data);

raw_data_avg_matrix = repmat(raw_data_avg, number_of_proteins, 1);

data_normalized = raw_data./raw_data_avg_matrix;

%% PCA 
% step 1: Calculate the average over all vectors
avg = mean(data_normalized);

%%
% step 2: Subtract the average from all vectors
diff_avg = data_normalized-repmat(avg, number_of_proteins, 1);

%%
% step 3: Calculate the covariance_matrixariance matrix

covariance_matrix = zeros(number_of_samples, number_of_samples);

for i = 1:number_of_samples;
    for j = 1:number_of_samples;
        covariance_matrix(i,j) = 1/(number_of_proteins - 1) * ...
            sum(diff_avg(:,i).*diff_avg(:,j));
    end
end
imagesc(covariance_matrix);
colorbar;
% disp(covariance_matrix);

%%
% step 4: calculate eigs of covariance_matrixariance matrix
[V,D] = eigs(covariance_matrix);
% disp(V);
% disp(diag(D));

%%
% plot data with the eigenvectors 
figure;
scatter3(data_normalized(:,1),data_normalized(:,2),...
    data_normalized(:,3),100,[1 0 0],'filled');
ev1 = V(:,1)*D(1,1);
ev2 = V(:,2)*D(2,2);
ev3 = V(:,3)*D(3,3);
hold on;
plot3([avg(1)-ev1(1) avg(1)+ev1(1)],[avg(2)-ev1(2) avg(2)+ev1(2)], ...
     [avg(3)-ev1(3) avg(3)+ev1(3)],'Color',[0 0 1],'LineWidth',3);
plot3([avg(1)-ev2(1) avg(1)+ev2(1)],[avg(2)-ev2(2) avg(2)+ev2(2)], ...
     [avg(3)-ev2(3) avg(3)+ev2(3)],'Color',[0 0 0],'LineWidth',3);
plot3([avg(1)-ev3(1) avg(1)+ev3(1)],[avg(2)-ev3(2) avg(2)+ev3(2)], ...
     [avg(3)-ev3(3) avg(3)+ev3(3)],'Color',[0 1 1],'LineWidth',3);
hold off;
axis equal

%%
 mapcaplot(data_normalized)






