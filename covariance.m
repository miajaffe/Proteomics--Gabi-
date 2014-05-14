%% Analyze covariance and correlation among replicates and among GI tract locations

clear all
close all hidden
clc
load('axes.mat');
load('OverlordMatrix.mat');
load('MusProt.mat')
normOverlord = OverlordNormalizer(2);
%%
% analyzing variance among mice 1,2,3 in the same colonization
% state/location along the GI tract. 
% The top 5 are from the GF mice, next 5 are B theta, and last 5 are
% conventional mice.
% The ordering from L to R: cecum, ileum, jejunum, prox colon, stomach 
figure;
count = 1;
for r = 1:3
    for s = 1:5
        cov_matrix = cov(normOverlord(:,:,r,s))
        subplot(3,5,count)
        count = count + 1;
        imagesc(cov_matrix)
    end
    suptitle('Variance Among Mice Replicates');
end

%%
% Correlation matrix for the same data
count = 1;
figure;
for r = 1:3
    for s = 1:5
        cor_matrix = corrcoef(normOverlord(:,:,r,s))
        subplot(3,5,count)
        count = count + 1;
        imagesc(cor_matrix)
    end
    suptitle('Correlation Among Mice Replicates');
end
%% Get mean and standard deviations of replicates - from Mia's scripts

% Calulate the mean and std deviation of each of the three replicates for
% each protein in each of the conditions, and store in two new 3-D
% matrices.

% loop through all 892 proteins
for i = 1:892
    % loop through all 3 colonization states
    for j = 1:3
      % loop through all 5 positions
      for k = 1:5
        % record mean and standard deviations to two new matrices
        AvgReps(i,j,k) = mean(normOverlord(i,:,j,k));
        StdReps(i,j,k) = std(normOverlord(i,:,j,k));
      end
   end
end

figure; 
% Analyze covariance of protein expression at different regions along the
% GI tract in GF Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,1,5) AvgReps(:,1,3) AvgReps(:,1,2), AvgReps(:,1,1), AvgReps(:,1,4)]
cov_matrix1 = cov(test)
subplot(3,1,1)
imagesc(cov_matrix1)
title('Covariance Matrix of samples along the GI tract, GF Mice')
% regions= {'Stomach'; 'Jejunum'; 'Ileum'; 'Cecum'; 'Proximal Colon'}
% xlabel(regions)

% Analyze covariance of protein expression at different regions along the
% GI tract in B theta colonized Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,2,5) AvgReps(:,2,3) AvgReps(:,2,2), AvgReps(:,2,1), AvgReps(:,2,4)]
cov_matrix2 = cov(test)
subplot(3,1,2)
imagesc(cov_matrix2)
title('Covariance Matrix of samples along the GI tract, B. theta Mice')

% Analyze covariance of protein expression at different regions along the
% GI tract in conventionally colonized Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,3,5) AvgReps(:,3,3) AvgReps(:,3,2), AvgReps(:,3,1), AvgReps(:,3,4)]
cov_matrix3 = cov(test)
subplot(3,1,3)
imagesc(cov_matrix3)
title('Covariance Matrix of samples along the GI tract, Conventional Mice')

%%
% Correlation
 
figure;
% Analyze correlation of protein expression at different regions along the
% GI tract in GF Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,1,5) AvgReps(:,1,3) AvgReps(:,1,2), AvgReps(:,1,1), AvgReps(:,1,4)]
cor_matrix1 =  corrcoef(test)
subplot(3,1,1)
imagesc(cor_matrix1)
title('Correlation Matrix of samples along the GI tract, GF Mice')
% regions= {'Stomach'; 'Jejunum'; 'Ileum'; 'Cecum'; 'Proximal Colon'}
% xlabel(regions)

% Analyze correlation of protein expression at different regions along the
% GI tract in B theta colonized Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,2,5) AvgReps(:,2,3) AvgReps(:,2,2), AvgReps(:,2,1), AvgReps(:,2,4)]
cor_matrix2 = corrcoef(test)
subplot(3,1,2)
imagesc(cor_matrix2)
title('Correlation Matrix of samples along the GI tract, B. theta Mice')

% Analyze correlation of protein expression at different regions along the
% GI tract in conventionally colonized Mice. All three replicates are averaged together. 
% The order is stomach, jejunum, ileum, cecum, proximal colon.
test = [AvgReps(:,3,5) AvgReps(:,3,3) AvgReps(:,3,2), AvgReps(:,3,1), AvgReps(:,3,4)]
cor_matrix3 = corrcoef(test)
subplot(3,1,3)
imagesc(cor_matrix3)
title('Correlation Matrix of samples along the GI tract, Conventional Mice')
