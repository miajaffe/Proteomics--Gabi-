%%IN PROGRESS

clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = OverlordNormalizer(2);

%analyzing variance among samples for cecum, GF
cov_matrix = cov(normOverlord(:,:,1,1))
imagesc(cov_matrix)

test = [normOverlord(:,1,1,1) normOverlord(:,1,1,2) normOverlord(:,1,1,3), normOverlord(:,1,1,4), normOverlord(:,1,1,5)]
cov_matrix2 = cov(test)
imagesc(cov_matrix2)





