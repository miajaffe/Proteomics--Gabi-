

% This code is essentially the same as the code for kmeans shown in
% lecture. I had to change the distance from sqEuclidean from correlation
% because the function was unable to cluster based on correlation.

% To generate figures, I ran this code 3 times once for all_matrices{1}=Germ-free, 
% once for all_matrices{2}=b. theta and once for all_matrices{3}=conventional.
cl = 8;
rows = cl/4;
sample = 2;

% Run kmeans analysis
 [cidx, ctrs] = kmeans(all_matrices{sample}, cl, 'dist', 'sqEuclidean', 'rep', 2000, 'disp', 'final');

% Plot clusters
figure
for c=1:cl
    subplot(rows,4,c);
    plot(1:5, all_matrices{sample}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering of 843 Host Secreted Proteins in B theta Mice')

bt_8_centroid = ctrs;
bt_8_clusters = cidx;


% 16 groups seems to be a good size, the clusters were identical the 2 times
% I ran the algorithm using the that many clusters. Using 1000 replicates
% also really improved the cluster identification.

%% Accessing protein IDs from cluster:
% cidx is a vector containing the assigned cluster for each protein.
% For example, if you wanted to pull out all proteins from cluster 5,
% simply use find(cidx==5) to get indexes of all proteins in that cluster.
% Then call the index i, in  axes{1,1}{1, i} to obtain the protein ID.
% or using 1 line of code: axes{1,1}{1,find(cidx==5)}

cl = 16;
rows = cl/4;
figure
for c=1:cl
    subplot(rows,4,c); plot(1:5, cv_16_centroid(c,:)'); axis tight
axis off
end
suptitle('Centroids')


%% Splitting proteins in each colonization state into high and low expression
% Will use data generated with 12 cluster for each state.
% Clusters 5 (56 proteins) and 7 (723 proteins) of Germ-Free are low
% Clusters 1 (671), 5 (82) and 10 (32) of B Theta are low
% Clusters 4 (680), 5 (32), 10 (16) and 11 (70) of conventional are low

% get indices of proteins with consistently low values and high values
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

%% Normalize high abundance proteins using zscore, and re-cluster
for i = 1:3
    z_high_matrices{i} = zscore(high_matrices{i}, 1, 2);
    subplot(2,2,i); plot(1:5, z_high_matrices{i})
end

%% Cluster using normalized abundances
%

cl = 7;
rows = cl/4;
sample = 1;

[cidx, ctrs] = kmeans(z_high_matrices{sample}, cl, 'dist', 'sqEuclidean', 'rep', 1000, 'disp', 'final');

% Plot clusters
figure
for c=1:cl
    subplot(3,3,c);
    plot(1:5, z_high_matrices{sample}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering')

high_gf_7_clusters = cidx;
high_gf_7_centroids = ctrs;

%% Pulling out protein IDS from most recent clusters
% To start with will pull out list of IDS from germ free and B theta where
% the expression is relatively high in the stomach compared to other
% locations.

% For germ-free this pattern is present in cluster 4, when clustering into
% 6 clusters
gf_ind_hi_stom = find(high_gf_6_clusters==4);


% For B theta this pattern was in cluster 3 when clustering into 6 clusters
bt_ind_hi_stom = find(high_bt_6_clusters==3);

% now need axes to reference which protein ids....
