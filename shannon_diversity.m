%% Shannon Diversity
clear all
close all hidden
clc
load('MusProt.mat');
load('axes140523.mat');
load('normOverlord_shannon.mat')
normOverlord = normOverlord_shannon * 1000;
%%
all_samples = [];
all_labels = {};
proteins = axes{1}'; %generates list of protein ids
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            all_samples = [all_samples normOverlord(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end
    end
end
%%
% mouse1 coloniz 1 loc 1-5
% mouse1 coloniz 2 loc 1-5
% mouse1 coloniz 3 loc 1-5
% mouse 2 coloniz 1 loc 1-5
% mouse 2 coloniz 2 loc 1-5
% mouse 2 coloniz 3 loc 1-5
% mouse 3 coloniz 1 loc 1-5
% mouse 3 coloniz 2 loc 1-5
% mouse 3 coloniz 3 loc 1-5

[H,VarH]=index_SaW(all_samples);
reshape1 = reshape(H, 5,9)'; %9 by 5 where each column is a different region along the GI tract
reshape1 = [reshape1(:,5) reshape1(:,3) reshape1(:,2) reshape1(:,1) reshape1(:,4)] %put in order of GI tract
coloniz_1  = mean(reshape1([1 4 7],:)); %takes average of replicates
sem_coloniz_1 = std(reshape1([1 4 7],:))./sqrt(3);

coloniz_2  = mean(reshape1([2 5 8],:));
sem_coloniz_2 = std(reshape1([2 5 8],:))./sqrt(3);

coloniz_3  = mean(reshape1([3 6 9],:));
sem_coloniz_3 = std(reshape1([3 6 9],:))./sqrt(3);
sem_matrix = [sem_coloniz_1' sem_coloniz_2' sem_coloniz_3'];
matrix_for_bar_graph = [coloniz_1' coloniz_2' coloniz_3']; %5 by 3 matrix where 5 is locations and 3 is the colonization states
%%
colormap('jet')
barweb(matrix_for_bar_graph, sem_matrix);
legend('Germ Free', 'B. Theta', 'Conventional')
title('Shannon Diversity of Samples')
ylabel('Shannon-Weiner Index')




    
