%% K-means clustering, eliminating low abundance clusters, and re-clustering high-abundance
%% proteins.
% By: Mia Jaffe
% This script has 3 parts:
% 1. It runs k-means analysis on the proteins contained in the all_matrices array generated
%    in the plot_abundance_v_location.m file. As part of this, it plots the clusters and 
%    centroids for a given run.
% 2. In the next section, low abundance proteins contained in clusters with a y-axis of 
%    less than 8e-3 were removed from the data set for each colonization state.
% 3. High abundance proteins were normalized with z-score and for each col. state were 
%    re-clustered again using k-means.
% 
% Data for the clusters is all saved and accessible in the mia_clustering_work_with_final_norm.mat
% file. Re-running the k-means will result re-arrangement in the cluster numbers, and most of
% this script will not work unless cluster numbers are changed within the script!

%% 1. K-means analysis
% This code is essentially the same as the code for kmeans shown in
% lecture. I had to change the distance from sqEuclidean from correlation
% because the function was unable to cluster based on correlation.

% To generate figures, I ran this code 9 times: trying 8, 12, and 16 clusters for each
% colonization state: all_matrices{1} = Germ-free, all_matrices{2} = b. theta and 
% all_matrices{3} = conventional.

cl = 8; % desired number of clusters
rows = cl/4; % assign number of rows for subplotting
sample = 2; % desired colonization state

% Run kmeans analysis with 2000 replicates
 [cidx, ctrs] = kmeans(all_matrices{sample}, cl, 'dist', 'sqEuclidean', 'rep', 2000, 'disp', 'final');

% Plot clusters
figure
for c=1:cl
    subplot(rows,4,c);
    plot(1:5, all_matrices{sample}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering of 843 Host Secreted Proteins in B theta Mice')

% assign variables based on which sample and how many clusters were run:
bt_8_centroid = ctrs;
bt_8_clusters = cidx;


%% A note on accessing protein IDs from cluster:
% cidx is a vector containing the assigned cluster for each protein.
% For example, if you wanted to pull out all proteins from cluster 5,
% simply use find(cidx==5) to get indexes of all proteins in that cluster.
% Then call the index i, in  axes{1,1}{1, i} to obtain the protein ID.
% or using 1 line of code: axes{1,1}{1,find(cidx==5)}

%% Plotting the centroids of the clusters
% I ran this chunk of code for any k-means analysis saved above which I wanted to visualize
% the centroids of the clusters. I had to adjust the cluster number, and the variable name
% in the plot function based on which run I wanted to plot.

cl = 16; % assign cluster number
rows = cl/4; % assign number of rows for subplot
figure
for c=1:cl
    subplot(rows,4,c); plot(1:5, cv_16_centroid(c,:)'); axis tight
axis off
end
suptitle('Centroids')

%% 2. Removing low abundance proteins from data set
% This section creates two matrices, one containing high abundance proteins
% and one containing low abundance proteins. The threshold was based on the 
% K-means analysis with 12 clusters, and "low abundance" clusters were defined as 
% having a y-axis less than 8e-3.

% The following clusters are below threshold, and must be separated:
% Clusters 5 (56 proteins) and 7 (723 proteins) of Germ-Free are low
% Clusters 1 (671), 5 (82) and 10 (32) of B Theta are low
% Clusters 4 (680), 5 (32), 10 (16) and 11 (70) of conventional are low

% get indices of proteins with consistently low values and high values from each col. state
i_gf = sort([find(gf_12_clusters==5); find(gf_12_clusters==7)]);
i_gf_hi = sort([find(gf_12_clusters>0 & gf_12_clusters<5); find(gf_12_clusters==6);find(gf_12_clusters>7 & gf_12_clusters<13)]);

i_bt = sort([find(bt_12_clusters==1); find(bt_12_clusters==5); find(bt_12_clusters==10)]);
i_bt_hi = sort([find(bt_12_clusters>1 & bt_12_clusters<5); find(bt_12_clusters>5 & bt_12_clusters<10); find(bt_12_clusters>10 & bt_12_clusters<13)]);

i_cv = sort([find(cv_12_clusters==4); find(cv_12_clusters==5); find(cv_12_clusters==10); find(cv_12_clusters==11)]);
i_cv_hi = sort([find(cv_12_clusters>0 & cv_12_clusters<4); find(cv_12_clusters>5 & cv_12_clusters<10); find(cv_12_clusters==12)]);

% construct new structure where only low proteins are included
low_matrices{1} = all_matrices{1}(i_gf,:);
low_matrices{2} = all_matrices{2}(i_bt,:);
low_matrices{3} = all_matrices{3}(i_cv,:);

% store high proteins in a matrix by removing rows in all_matrices containing
% low abundance proteins
high_matrices{1} = removerows(all_matrices{1}, 'ind', i_gf);
high_matrices{2} = removerows(all_matrices{2}, 'ind', i_bt);
high_matrices{3} = removerows(all_matrices{3}, 'ind', i_cv);

%% 3. Normalize high abundance proteins using zscore, and re-cluster
for i = 1:3
    z_high_matrices{i} = zscore(high_matrices{i}, 1, 2);
    subplot(2,2,i); plot(1:5, z_high_matrices{i})
end

%% Cluster using normalized abundances
% This section was also run 9 times, once for cluster sizes of 5, 6 and 7 using data for 
% each of the three colonization states.

cl = 7; % set cluster number
rows = cl/4;
sample = 1; % set colonization state

% run k-means
[cidx, ctrs] = kmeans(z_high_matrices{sample}, cl, 'dist', 'sqEuclidean', 'rep', 1000, 'disp', 'final');

% Plot clusters
figure
for c=1:cl
    subplot(3,3,c);
    plot(1:5, z_high_matrices{sample}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering')

% save cluster information in chosen variable names
high_gf_7_clusters = cidx;
high_gf_7_centroids = ctrs;
