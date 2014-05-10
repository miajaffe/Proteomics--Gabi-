

% This code is essentially the same as the code for kmeans shown in
% lecture. I had to change the distance from sqEuclidean from correlation
% because the function was unable to cluster based on correlation.

% To generate figures, I ran this code 3 times once for all_matrices{1}=Germ-free, 
% once for all_matrices{2}=b. theta and once for all_matrices{3}=conventional.

% Run kmeans analysis
[cidx, ctrs] = kmeans(all_matrices{3}, 9, 'dist', 'sqEuclidean', 'rep', 5, 'disp', 'final');

% Plot clusters
figure
for c=1:9
    subplot(3,3,c);
    plot(1:5, all_matrices{3}((cidx ==c),:)'); axis tight
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'Stomach', 'Jejunum', 'Ileum', 'Cecum', 'Colon'})
end
suptitle('K-Means Clustering of 892 Host Secreted Proteins in Conventional Mice')