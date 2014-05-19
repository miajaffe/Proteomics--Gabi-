clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = OverlordNormalizer(2); %normalize by row
normOverlord = normOverlord .* 100;

cecum_samples = [];
for j = 1:3    
    samples = normOverlord(:,:,j,1)';
    cecum_samples = [cecum_samples; samples];
end
var_cecum_samples = var(cecum_samples);

cecum_GF = normOverlord(:,:,1,1)';
var_cecum_GF = var(cecum_GF);
var_cecum_BT = var(normOverlord(:,:,2,1)');
var_cecum_conv = var(normOverlord(:,:,3,1)');


    