% Examines incidence of proteins and/or GO codes of interest within
% normOverlordFinal and GOenrichMat
clear all
close all
% clc
load GOenrichMat
load GOenrichMat_shannon
load axes
load GOtoIndexConverterStr
load IndextoGOConverterStr
load allGODic
% %% Examine for REG3G O09049
% %http://stackoverflow.com/questions/8061344/how-to-search-for-a-string-in-cell-array-in-matlab
% protInd = find(ismember(axes{1},'O09049'));
% GF_BT_RF_prot_prevalence(1) = sum(sum(GOenrichMat_shannon(protInd,:,1,:)));
% GF_BT_RF_prot_prevalence(2) = sum(sum(GOenrichMat_shannon(protInd,:,2,:)));
% GF_BT_RF_prot_prevalence(3) = sum(sum(GOenrichMat_shannon(protInd,:,3,:)));
% %% Examine for GO:0050830: defense response to Gram-positive bacterium
% GO = 'GO:0050830';
% for i = 1:1:3
%     index = GOtoIndexConverterStr(GO);
%     GF_BT_RF_GO1(1,i) = sum(GOenrichMat_shannon(index,i,1,:),4);
%     GF_BT_RF_GO1(2,i) = sum(GOenrichMat_shannon(index,i,2,:),4);
%     GF_BT_RF_GO1(3,i) = sum(GOenrichMat_shannon(index,i,3,:),4);
% end
% GF_BT_RF_GO1mean = mean(GF_BT_RF_GO1,2);
% % t-tests
% % GF vs BT
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(2,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and BT with p-value: %f\n',GO,p);
% % BT vs RF
% [h,p] = ttest2(GF_BT_RF_GO1(2,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between BT and RF with p-value: %f\n',GO,p);
% % GF and RF
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and RF with p-value: %f\n',GO,p);
% %% Examine for GO:0050829: defense response to Gram-negative bacterium
% GO = 'GO:0050829';
% for i = 1:1:3
%     index = GOtoIndexConverterStr(GO);
%     GF_BT_RF_GO1(1,i) = sum(GOenrichMat_shannon(index,i,1,:),4);
%     GF_BT_RF_GO1(2,i) = sum(GOenrichMat_shannon(index,i,2,:),4);
%     GF_BT_RF_GO1(3,i) = sum(GOenrichMat_shannon(index,i,3,:),4);
% end
% GF_BT_RF_GO1mean = mean(GF_BT_RF_GO1,2);
% % t-tests
% % GF vs BT
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(2,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and BT with p-value: %f\n',GO,p);
% % BT vs RF
% [h,p] = ttest2(GF_BT_RF_GO1(2,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between BT and RF with p-value: %f\n',GO,p);
% % GF and RF
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and RF with p-value: %f\n',GO,p);
% 
% %% Examine for GO:0000272: polysaccharide catabolic process
% GO = 'GO:0000272';
% for i = 1:1:3
%     index = GOtoIndexConverterStr(GO);
%     GF_BT_RF_GO1(1,i) = sum(GOenrichMat_shannon(index,i,1,:),4);
%     GF_BT_RF_GO1(2,i) = sum(GOenrichMat_shannon(index,i,2,:),4);
%     GF_BT_RF_GO1(3,i) = sum(GOenrichMat_shannon(index,i,3,:),4);
% end
% GF_BT_RF_GO1mean = mean(GF_BT_RF_GO1,2);
% % t-tests
% % GF vs BT
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(2,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and BT with p-value: %f\n',GO,p);
% % BT vs RF
% [h,p] = ttest2(GF_BT_RF_GO1(2,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between BT and RF with p-value: %f\n',GO,p);
% % GF and RF
% [h,p] = ttest2(GF_BT_RF_GO1(1,:),GF_BT_RF_GO1(3,:),'Vartype','unequal');
% fprintf('Incidence of %s is different between GF and RF with p-value: %f\n',GO,p);
%% Find GO codes with p-value <= 0.01/2991
allGO = keys(GOtoIndexConverterStr);
for i = 1:1:length(allGO)
    GOcurr = allGO{i};
    for j = 1:1:3
        index = GOtoIndexConverterStr(GOcurr);
        GF_BT_RF_GO1(1,j) = sum(GOenrichMat_shannon(index,j,1,:),4);
        GF_BT_RF_GO1(2,j) = sum(GOenrichMat_shannon(index,j,2,:),4);
        GF_BT_RF_GO1(3,j) = sum(GOenrichMat_shannon(index,j,3,:),4);
    end
    p1 = anova1(horzcat(GF_BT_RF_GO1(1,:)',GF_BT_RF_GO1(2,:)'),{'GF','BT'},'off');
    p2 = anova1(horzcat(GF_BT_RF_GO1(2,:)',GF_BT_RF_GO1(3,:)'),{'BT','RF'},'off');
    p3 = anova1(horzcat(GF_BT_RF_GO1(1,:)',GF_BT_RF_GO1(3,:)'),{'GF','RF'},'off');
    p(i,1) = p1;
    p(i,2) = p2;
    p(i,3) = p3;
    pmin(i) = min(p(i,:));
end
[sortedpmin, pindeces] = sort(pmin,'ascend');
sortedp = p(pindeces,:);
allGOsorted = allGO(pindeces);
counter = 1;
while sortedpmin(counter) < 0.01/2991
    value = allGOsorted(counter);
    temp = allGODic(value{1});
    significantGO{counter,1} = temp{1};
    significantGO{counter,2} = sortedp(counter,1);
    significantGO{counter,3} = sortedp(counter,2);
    significantGO{counter,4} = sortedp(counter,3);
    counter = counter + 1;
end
fileID = fopen('significantGOtitles.dat','w');
formatSpec0 = '%s \t %s \t %s \t %s\n\n\n';
header = {'GO Title','GF-BT diff','BT-RF diff','GF-RF diff'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%s\t%1.4f\t%1.4f\t%1.4f\n';
[nrows, ncols] = size(significantGO);
for row = 1:nrows
    fprintf(fileID,formatSpec,significantGO{row,:});
end
fclose(fileID);
% Look for GO terms for which RF is significantly different from both GF
% and BT
counter = 1;
significantGORF = {};
for i = 1:1:size(sortedp,1)
    if sortedp(i,2) <= 0.01/2991 && sortedp(i,3) <= 0.01/2991
        value = allGOsorted(i);
        temp = allGODic(value{1});
        significantGORF{counter,1} = temp{1};
        significantGORF{counter,2} = sortedp(i,1);
        significantGORF{counter,3} = sortedp(i,2);
        significantGORF{counter,4} = sortedp(i,3);
        counter = counter + 1;
    end
end
fileID = fopen('GOdiffFromRF.dat','w');
formatSpec0 = '%s \t %s \t %s \t %s\n\n\n';
header = {'GO Title','GF-BT diff','BT-RF diff','GF-RF diff'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%s\t%1.4f\t%1.4f\t%1.4f\n';
[nrows, ncols] = size(significantGORF);
for row = 1:nrows
    fprintf(fileID,formatSpec,significantGORF{row,:});
end
fclose(fileID);
    