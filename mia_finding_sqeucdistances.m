% CALCULATE mean squared euclidean distance between centroids


% for all 6 CV centroids
for i = 1:6
    % for all 5 GF centroids
    for j = 1:5
        % find mean squared eucledian distance
        c = high_cv_6_centroids(i,:);
        b = high_gf_5_centroids(j,:);
        msed = mean(abs(b-c).^2);
        msed_mat(i,j) = msed
        % store in matrix
    end
end
