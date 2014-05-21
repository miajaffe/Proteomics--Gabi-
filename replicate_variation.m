clear all
close all hidden
clc
load('axes.mat');
load('LetterMap.mat');
load('normOverlord2.mat');
load('ProteinMap.mat');
load('MusProt.mat')
normOverlord = normOverlord2 .* 100;
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
hold on;
x = linspace(0,0.25,100);
y = x;
plot(x,y, 'r')


%replicate 1 vs 3
figure
scatter(log(replicate_one), log(replicate_three))
xlabel('Replicate 1')
ylabel('Replicate 3')


%replicate 2 vs 3
figure 
scatter(log(replicate_two), log(replicate_three))
xlabel('Replicate 2')
ylabel('Replicate 3')

%%
%3D plot of all variables
figure
scatter3(log(replicate_one), log(replicate_two),log(replicate_three));



        
