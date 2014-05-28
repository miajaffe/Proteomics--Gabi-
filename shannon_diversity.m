%% Shannon Diversity
clear all
close all hidden
clc
load('MusProt.mat');
load('axes140523.mat');
load('normOverlordFinal_140523.mat');
load('OverlordMatrix');

normOverlord_manual = OverlordNormalizer();
normOverlord = normOverlordFinal * 1000;
%%
cutoff = 5;
cutoff_matrix = OverlordMatrix;
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            [rows] = find(OverlordMatrix <= cutoff);
            cutoff_matrix(rows,1) = 0;
        end
    end
end
%%
%[H,VarH]=index_SaW(OverlordMatrix,2);

[H,VarH]=index_SaW(normOverlord,2);
    
