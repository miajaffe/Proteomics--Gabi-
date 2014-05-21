

% This code is essentially the same as the code for kmeans shown in
% lecture. I had to change the distance from sqEuclidean from correlation
% because the function was unable to cluster based on correlation.

% To generate figures, I ran this code 3 times once for all_matrices{1}=Germ-free, 
% once for all_matrices{2}=b. theta and once for all_matrices{3}=conventional.
cl = 8;
rows = cl/4;
sample = 3;

% Run kmeans analysis
 [cidx, ctrs] = kmeans(all_matrices{sample}, cl, 'dist', 'sqEuclidean', 'rep', 2000, 'disp', 'final');

% Plot clusters
figure
for c=1:cl
    subplot(rows,4,c);
    plot(1:5, all_matrices{sample}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'St', 'Jej', 'Ile', 'Cec', 'Col'})
end
suptitle('K-Means Clustering of 876 Host Secreted Proteins in Conventional Mice')

cv_8_centroid = ctrs;
cv_8_clusters = cidx;


% 16 groups seems to be a good size, the clusters were identical the 2 times
% I ran the algorithm using the that many clusters. Using 1000 replicates
% also really improved the cluster identification.

%% Accessing protein IDs from cluster:
% cidx is a vector containing the assigned cluster for each protein.
% For example, if you wanted to pull out all proteins from cluster 5,
% simply use find(cidx==5) to get indexes of all proteins in that cluster.
% Then call the index i, in  axes{1,1}{1, i} to obtain the protein ID.
% or using 1 line of code: axes{1,1}{1,find(cidx==5)}

cl = 8;
rows = cl/4;
figure
for c=1:cl
    subplot(rows,4,c); plot(1:5, cv_8_centroid(c,:)'); axis tight
axis off
end
suptitle('Centroids')

