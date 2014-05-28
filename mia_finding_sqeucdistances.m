%% Identifying distances between clusters and proteins within clusters
% By: Mia Jaffe
% To run this script, make sure to load relevant data in mia_clustering_work_with_final_norm.mat
% This script has several parts:
% 1. It calculates the squared euclidean distance between clusters generated using the
%    high abundance proteins generated in k_means.
% 2. This section allows cluster numbers for each of the 3 states to entered. Then, it
%    retrieves and prints the number of and IDs of the proteins in the given clusters. It
%    also plots the original traces of the proteins across locations prior to z-score
%    normalization in order to verify the correct cluster is being reported.
% 3. This section is run when the IDs and number of overlapping proteins is desired between
%    the GF and BT clusters entered in section 2.
% 4. This section is run when the IDs and number of overlapping proteins is desired between
%    the GF and CV clusters entered in section 2.
% 5. This section was used to determine the total overlap in protein IDs between two 
%    colonization states.


%% 1. Calculate mean squared euclidean distance between centroids
% This section was run twice such that the centroids of each of the 6 conventional 
% clusters, and each of the 5 B. theta clusters were compared to the centroids of each of 
% the 5 Germ-Free clusters. It was done this way because we decided to use Germ-free as a
% baseline of analysis. Results from the command window were copied to separate data files.

% for all 5 GF centroids
for i = 1:5
    % for all 6 BT centroids
    for j = 1:6
        % find mean squared eucledian distance
        c = high_gf_5_centroids(i,:);
        b = high_cv_6_centroids(j,:);
        msed = sum(abs(b-c).^2);
        msed_mat(i,j) = msed
        % store in matrix
    end
end


%% Find the number of and names of proteins in desired cluster:

% Use 'high_cv_6_clusters' for conventional
% Use 'high_bt_5_clusters' for b theta
% Use 'high_gf_5_clusters' for germ-free

% 'i_bt_hi' has normOverlord indices of reclustered proteins 
% 'i_gf_hi' and 'i_cv_hi' have the same for germ-free and conventional

% enter the cluster number from which to retrieve protein information

desired_cluster_bt = 1;
desired_cluster_gf = 1;
desired_cluster_cv = 1;

% Get B. theta proteins ids in cluster
cluster_bt = find(high_bt_5_clusters==desired_cluster_bt);
original_i_bt = i_bt_hi(cluster_bt);
ID_bt = axes{1,1}(1, original_i_bt);


% Get Germ-free proteins ids in cluster
cluster_gf = find(high_gf_5_clusters==desired_cluster_gf);
original_i_gf = i_gf_hi(cluster_gf);
ID_gf = axes{1,1}(1, original_i_gf);

% Get conventional proteins ids in cluster:
cluster_cv = find(high_cv_6_clusters==desired_cluster_cv);
original_i_cv = i_cv_hi(cluster_cv);
ID_cv = axes{1,1}(1, original_i_cv);


% plot original values to determine if correct traces were pulled from
% clusters:
subplot(2,2,1); plot(1:5, (all_matrices{2}((original_i_bt),:)'))
subplot(2,2,2); plot(1:5, (all_matrices{1}((original_i_gf),:)'))
subplot(2,2,3); plot(1:5, (all_matrices{3}((original_i_cv),:)'))

% Print out number of and IDs for proteins in cluster for each colonization
% state:


disp('Number of B theta proteins: ')
disp(length(ID_bt))

disp('IDs of B theta proteins: ')
for i = 1:length(ID_bt)
   disp(ID_bt{i})
end


disp('Number of germ-free proteins: ')
disp(length(ID_gf))

disp('IDs of germ-free proteins: ')
for i = 1:length(ID_gf)
    disp(ID_gf{i})
end

disp('Number of conventional proteins: ')
disp(length(ID_cv))

disp('IDs of conventional proteins: ')
for i = 1:length(ID_cv)
    disp(ID_cv{i})
end

%% 3. Compare IDs within the clusters assigned above in B theta vs Germ-Free
% Use this section when the overlap between BT and GF clusters is desired.

disp('Proteins found in both germ-free and B theta: ')

count = 0;
for k = 1:length(ID_bt)
   for m = 1:length(ID_gf) 
      gf_id = ID_gf{m};
      bt_id = ID_bt{k};
       c = strcmp(ID_bt{k}, ID_gf{m});
       if c == 1
           disp(gf_id)
           
           %fprintf('%s in B theta matches %s in Germ Free\n', bt_id, gf_id)
           count = count +1;
       end
   end
end

disp(count)

%% 4. Compare IDs within the clusters assigned above in Conventional vs Germ-Free
% Use this section when the overlap between CV and GF clusters is desired.

disp('Proteins found in both germ-free and conventional: ')

count = 0;
for m = 1:length(ID_cv)
   for n = 1:length(ID_gf) 
      gf_id = ID_gf{n};
      cv_id = ID_cv{m};
       c = strcmp(ID_cv{m}, ID_gf{n});
       if c == 1
           disp(gf_id)
           %fprintf('%s in B theta matches %s in Germ Free\n', bt_id, gf_id)
           count = count +1;
       end
   end
end

disp(count)

%% Compare all IDs contained in filtered data between germ-free and conventional
count = 0;

for i = 1:length(i_gf_hi)
    for j = 1:length(i_cv_hi)
        if i_gf_hi(i) == i_cv_hi(j)
            count = count+1;
        end
        
    end
end
disp(count)
