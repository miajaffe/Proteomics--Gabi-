clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('OverlordMatrix.mat');
load('normOverlordFinal.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = normOverlordFinal .* 100;
%%
replicate_one = [];
replicate_two = [];
replicate_three = [];
for colonization = 1:3
    for loc = 1:5
        replicate_one = [replicate_one; normOverlord(:,1, colonization, loc)];
        replicate_two = [replicate_two; normOverlord(:,2, colonization, loc)];
        replicate_three = [replicate_three; normOverlord(:,3, colonization, loc)];
    end
end
%%
%replicate 1 vs 2
scatter(log(replicate_one), log(replicate_two))
xlabel('Replicate 1')
ylabel('Replicate 2')
%plot y = x
%hold on;
%x = linspace(0,0.25,100);
%y = x;
%plot(x,y, 'r')
%Get rid of proteins with 0 abundance in either replicate
indices_zeros_rep1 = find(replicate_one == 0);
indices_zeros_rep2 = find(replicate_two == 0);
indices_zeros_rep3 = find(replicate_three == 0);
zeros_total_indices = [indices_zeros_rep1; indices_zeros_rep2];
zeros_total_indices = unique(zeros_total_indices);
rep_one_removed = replicate_one;
rep_one_removed(zeros_total_indices) = [];
rep_two_removed = replicate_two;
rep_two_removed(zeros_total_indices) = [];

reg_1 = fitlm(log(rep_one_removed),log(rep_two_removed))


%replicate 1 vs 3
figure
scatter(log(replicate_one), log(replicate_three));
xlabel('Replicate 1')
ylabel('Replicate 3')
zeros_total_indices_2 = [indices_zeros_rep1; indices_zeros_rep3];
zeros_total_indices_2 = unique(zeros_total_indices_2);
rep_one_removed = replicate_one;
rep_one_removed(zeros_total_indices_2) = [];
rep_three_removed = replicate_three;
rep_three_removed(zeros_total_indices_2) = [];
reg_2 = fitlm(log(rep_one_removed),log(rep_three_removed))


%replicate 2 vs 3
figure 
scatter(log(replicate_two), log(replicate_three))
xlabel('Replicate 2')
ylabel('Replicate 3')
zeros_total_indices_3 = [indices_zeros_rep2; indices_zeros_rep3];
zeros_total_indices_3 = unique(zeros_total_indices_3);
rep_two_removed = replicate_two;
rep_two_removed(zeros_total_indices_3) = [];
rep_three_removed = replicate_three;
rep_three_removed(zeros_total_indices_3) = [];
reg_3 = fitlm(log(rep_two_removed),log(rep_three_removed))

%%
%3D plot of all variables
figure
scatter3(log(replicate_one), log(replicate_two),log(replicate_three));



        
