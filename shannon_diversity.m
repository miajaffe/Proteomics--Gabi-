%% Shannon Diversity
clear all
close all hidden
clc
load('axes140523.mat');
load('normOverlord_shannon.mat')
load('normOverlordFinal_140523.mat');

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

%% t-tests
GF_cecum = reshape1([1 4 7],4);
GF_prox_colon= reshape1([1 4 7],5);
[H, P] = ttest(GF_cecum, GF_prox_colon);
[P_rank, H_rank] = ranksum(GF_cecum, GF_prox_colon);

BT_cecum = reshape1([2 5 8],4);
BT_prox_colon= reshape1([2 5 8],5);
[H2, P2] = ttest(BT_cecum, BT_prox_colon);
[P2_rank, H2_rank] = ranksum(BT_cecum, BT_prox_colon);

%% anova - cecum
anova_matrix_cecum = [reshape1([1 4 7],4) reshape1([2 5 8],4) reshape1([3 6 9],4)];
[p_anova, table, stats] = anova1(anova_matrix_cecum);
figure
sig_cecum = multcompare(stats)
%% anova - prox colon
anova_matrix_colon = [reshape1([1 4 7],5) reshape1([2 5 8],5) reshape1([3 6 9],5)];
[p_anova2, table2, stats2] = anova1(anova_matrix_colon);
figure
sig_colon = multcompare(stats2)

%%
colormap('jet')
subplot(2,1,1)
bar_graph = barweb(matrix_for_bar_graph, sem_matrix);
% set(bar_graph,'FaceColor','r');
legend('Germ Free', 'B. Theta', 'Conventional')
title('Shannon Diversity of Samples')
ylabel('Shannon-Weiner Index')


%%

%%calculate the number of unique proteins in each sample
protein_num= zeros(45,1);
for sample = 1:45
   non_zero_indices = find(all_samples(:,sample));
   prot_num(sample) = length(non_zero_indices);
end
prot_num_reshaped = reshape(prot_num, 5,9)';
prot_num_reshaped = [prot_num_reshaped(:,5) prot_num_reshaped(:,3) prot_num_reshaped(:,2) prot_num_reshaped(:,1) prot_num_reshaped(:,4)]; %put in order of GI tract
coloniz_1_prot_num  = mean(prot_num_reshaped([1 4 7],:)); %takes average of replicates
sem_coloniz_1_prot = std(prot_num_reshaped([1 4 7],:))./sqrt(3);
coloniz_2_prot_num  = mean(prot_num_reshaped([2 5 8],:)); %takes average of replicates
sem_coloniz_2_prot = std(prot_num_reshaped([2 5 8],:))./sqrt(3);
coloniz_3_prot_num  = mean(prot_num_reshaped([3 6 9],:)); %takes average of replicates
sem_coloniz_3_prot = std(prot_num_reshaped([3 6 9],:))./sqrt(3);

sem_matrix_prot = [sem_coloniz_1_prot' sem_coloniz_2_prot' sem_coloniz_3_prot'];
protein_num_matrix = [coloniz_1_prot_num' coloniz_2_prot_num' coloniz_3_prot_num'];%5 by 3 matrix where 5 is locations and 3 is the colonization states

subplot(2,1,2)
bar_graph = barweb(protein_num_matrix, sem_matrix_prot);
legend('Germ Free', 'B. Theta', 'Conventional')
ylabel('Unique Proteins Identified') 

%% anova - cecum - number proteins - GF has significantly more proteins than conventional
anova_cecum_nprot = [prot_num_reshaped([1 4 7],4) prot_num_reshaped([2 5 8],4) prot_num_reshaped([3 6 9],4)];
[p_anova3, table3, stats3] = anova1(anova_cecum_nprot);
figure
sig_cecum_nprot = multcompare(stats3)
%% anova - prox colon - number proteins - so so close to GF having significantly more proteins
anova_colon_nprot = [prot_num_reshaped([1 4 7],5) prot_num_reshaped([2 5 8],5) prot_num_reshaped([3 6 9],5)];
[p_anova4, table4, stats4] = anova1(anova_colon_nprot);
figure
sig_colon_nprot = multcompare(stats4)

    
