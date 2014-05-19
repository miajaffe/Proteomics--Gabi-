

% This code is essentially the same as the code for kmeans shown in
% lecture. I had to change the distance from sqEuclidean from correlation
% because the function was unable to cluster based on correlation.

% To generate figures, I ran this code 3 times once for all_matrices{1}=Germ-free, 
% once for all_matrices{2}=b. theta and once for all_matrices{3}=conventional.

% Run kmeans analysis
 [cidx, ctrs] = kmeans(norm_col_all_matrices{2}, 16, 'dist', 'sqEuclidean', 'rep', 1000, 'disp', 'final');

% Plot clusters
figure
for c=1:16
    subplot(4,4,c);
    plot(1:5, norm_col_all_matrices{2}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering of 892 Host Secreted Proteins in Germ-Free Mice')

% 16 groups seems to be a good size, the clusters were identical the 2 times
% I ran the algorithm using the that many clusters. Using 1000 replicates
% also really improved the cluster identification.

%% Accessing protein IDs from cluster:
% cidx is a vector containing the assigned cluster for each protein.
% For example, if you wanted to pull out all proteins from cluster 5,
% simply use find(cidx==5) to get indexes of all proteins in that cluster.
% Then call the index i, in  axes{1,1}{1, i} to obtain the protein ID.
% or using 1 line of code: axes{1,1}{1,find(cidx==5)}


figure
for c=1:16
    subplot(4,4,c); plot(1:5,ctrs3(c,:)'); axis tight
axis off
end
suptitle('K-Means Clustering of 310 Yeast Genes, Centroids')

