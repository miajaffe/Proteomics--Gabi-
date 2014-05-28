%% Data preparation and plotting protein abundance by location
% By: Mia Jaffe
% This script has 4 parts:
% 1. It assigns which matrix to pull data from
% 2. It calculates the mean and standard deviation of the 3 mouse replicates
% 3. It plots, for each protein in each colonization state the abundance vs location. 
%    This section of the code also generates a new array, all_matrices, which contains
%    one matrix for each colonization state holding avg counts of proteins (rows) 
%    at each location (columns). This matrix will be referenced in the k-means script.
% 4. Finally, it has the ability to plot the proteins vs locations excluding proteins 
%    above and below a given threshold. This section was used in initial tests. 

%% Call normalization script to normalize values in matrix

% NewOverlord = normalization_mia(OverlordMatrix);
% 5/20/14-no longer needed, will use Evan's normalized data in
% normOverlordF.

NewOverlord = normOverlordFinal;

%% Get mean and standard deviations of replicates

% Calulate the mean and std deviation of each of the three replicates for
% each protein in each of the conditions, and store in two new 3-D
% matrices.

% loop through all proteins
for i = 1:843
    % loop through all 3 colonization states
    for j = 1:3
      % loop through all 5 positions
      for k = 1:5
        % record mean and standard deviations to two new matrices
        AvgReps(i,j,k) = mean(NewOverlord(i,:,j,k));
        StdReps(i,j,k) = std(NewOverlord(i,:,j,k));
      end
   end
end




%% Plot line graph with protein abundance (y-axis) vs. location (x-axis)

% Currently indices of gut locations are in the following order:
% 1 = cecum
% 2 = ileum
% 3 = jejunem
% 4 = proximal colon
% 5 = stomach

% Would like to plot in the following order:
% 1 = stomach
% 2 = jejunum
% 3 = ileum
% 4 = cecum
% 5 = proximal colon

% initiate variable to store colonization state names (to be used for plot
% titles)
colStateKey = {'Germ-Free', 'B. theta', 'Conventional'};

clf;

% loop through three colonization states:
for i = 1:3
    
    % record protein values in each location
    protCec = AvgReps(:,i,1); 
    protIle = AvgReps(:,i,2);
    protJej = AvgReps(:,i,3);
    protCol = AvgReps(:,i,4);
    protStom = AvgReps(:,i,5);
    
    % make new matrix where rows are proteins, and columns are locations
    newM = [protStom, protJej, protIle, protCec,protCol]; 
    
    % Plot abundance over location
    subplot(2,2,i); plot(newM')
    
    % add title with colonization state stated
    title(sprintf('Abundance Across Locations for %s', colStateKey{i}));
    ylabel('Count')
    
    % change x-axis labels to be the locations
    set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'Stomach', 'Jejunum', 'Ileum', 'Cecum', 'Colon'})
     
    % save matrix generated in this interation to a cell array that will 
    % contain the matrices for all three colonization states:
    all_matrices{i} = newM;
      
end

%% For easier viewing, re-plot excluding proteins with very high or very low abundance

% For each matrix (one for each colonization state), calculate the mean of 
% each protein over all 5 locations, and make a new matrix containing only 
% proteins with a mean abundance between the set thresholds.

% set minimum and maximum protein abundance averages:
max = .100;
min = .010;

clf;

% loop through 3 colonization states
for i = 1:3
    
   % find avg count for each protein over 5 locations/columns:
   row_means = mean(all_matrices{i},2);
   
   % make a new matrix of proteins with average expression
   % greater than the threshold:
   matrix_above_thresh = all_matrices{i}(row_means > min, :);
   
   % find row means again for new matrix:
   row_means2 = mean(matrix_above_thresh,2);
   
   % make new matrix of proteins with average expression below threshold:
   matrix_between_thresh = matrix_above_thresh((row_means2 < max), :);
   
   % plot!
   subplot(2,2,i); plot(matrix_between_thresh')
   ylim([0 200])
   
   % add title with colonization state stated
   title(sprintf('Abundance Across Locations for %s', colStateKey{i}));
   ylabel('Count')
    
   % change x-axis labels to be the locations
   set(gca,'XTick', [1 2 3 4 5], 'XTickLabel', {'Stomach', 'Jejunum', 'Ileum', 'Cecum', 'Colon'})
    
    
end









