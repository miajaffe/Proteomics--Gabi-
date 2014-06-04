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
load titleToGO
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
%% Find GO codes with q-value <= arbitrary such that less than 1 significant sample is a false positive
allGO = keys(GOtoIndexConverterStr);
for i = 1:1:length(allGO)
    GOcurr = allGO{i};
    for j = 1:1:3
        index = GOtoIndexConverterStr(GOcurr);
        GF_BT_RF_GO1(j,1) = sum(GOenrichMat_shannon(index,j,1,:),4);
        GF_BT_RF_GO1(j,2) = sum(GOenrichMat_shannon(index,j,2,:),4);
        GF_BT_RF_GO1(j,3) = sum(GOenrichMat_shannon(index,j,3,:),4);
    end
    [p,table,stats] = anova1(GF_BT_RF_GO1,{'GF','BT','RF'},'off');
%     c(i,:,:) = multcompare(stats);
    pall(i) = p;
    statsall(i) = stats;
end
pall_sorted = sort(pall,'ascend');
[FDR, q] = mafdr(pall);
[sortedq, qind] = sort(q,'ascend');
counter = 1;
qcutoff = 0.05;
statssorted = statsall(qind);
while sortedq(counter) <= qcutoff
    tempnum = mod(qind(counter),3);
    index = qind(counter);
    value = allGO(index);
    temp = allGODic(value{1});
    counts = titleToGO(temp{1});
    counts = counts{3};
    significantGO{counter,1} = temp{1};
    significantGO{counter,2} = mean(sum(counts(1,:,1,:),4));
    significantGO{counter,3} = mean(sum(counts(1,:,2,:),4));
    significantGO{counter,4} = mean(sum(counts(1,:,3,:),4));
    significantGO{counter,5} = std(sum(counts(1,:,1,:),4));
    significantGO{counter,6} = std(sum(counts(1,:,2,:),4));
    significantGO{counter,7} = std(sum(counts(1,:,3,:),4));
    significantGO{counter,8} = pall_sorted(counter);
    significantGO{counter,9} = sortedq(counter);
    allmeans(counter) = mean([mean(sum(counts(1,:,1,:),4)),mean(sum(counts(1,:,2,:),4)),mean(sum(counts(1,:,3,:),4))]);
    counter = counter + 1;
end
[sortedallmeans, meanind] = sort(allmeans,'descend');
statsfinal = statssorted(meanind);
fileID = fopen('significantGOtitlesQvalsCorrected5pctFDR.csv','w');
formatSpec0 = '%s , %s , %s , %s, %s, %s, %s, %s, %s\n';
header = {'GO Title','GF count','BT count','RF count','GF std','BT std','RF std','p-val','q-val'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%s,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n';
[nrows, ncols] = size(significantGO);
for row = 1:nrows
    fprintf(fileID,formatSpec,significantGO{row,:});
end
fclose(fileID);
for j = 1:1:15
    figure
    multcompare(statsfinal(j),'alpha',0.01);
end
% %% Same analysis as above but for GI location
% for i = 1:1:length(allGO)
%     GOcurr = allGO{i};
%     for j = 1:1:3
%         index = GOtoIndexConverterStr(GOcurr);
%         GUT_GO1(j,1) = sum(GOenrichMat_shannon(index,j,:,1),3);
%         GUT_GO1(j,2) = sum(GOenrichMat_shannon(index,j,:,2),3);
%         GUT_GO1(j,3) = sum(GOenrichMat_shannon(index,j,:,3),3);
%         GUT_GO1(j,4) = sum(GOenrichMat_shannon(index,j,:,4),3);
%         GUT_GO1(j,5) = sum(GOenrichMat_shannon(index,j,:,5),3);
%     end
%     [p,table,stats] = anova1(GUT_GO1,{'cecum','ilum','jejunum','prox colon','Stomach'},'off');
%     c(i) = multcompare(stats);
%     pall(i) = p;
% end
% pall_sorted = sort(pall,'ascend');
% [FDR, q] = mafdr(pall);
% [sortedq, qind] = sort(q,'ascend');
% counter = 1;
% qcutoff = 0.003;
% while sortedq(counter) <= qcutoff
%     tempnum = mod(qind(counter),3);
%     index = qind(counter);
%     value = allGO(index);
%     temp = allGODic(value{1});
%     counts = titleToGO(temp{1});
%     counts = counts{3};
%     significantGO{counter,1} = temp{1};
%     significantGO{counter,2} = mean(sum(counts(1,:,:,1),3));
%     significantGO{counter,3} = mean(sum(counts(1,:,:,2),3));
%     significantGO{counter,4} = mean(sum(counts(1,:,:,3),3));
%     significantGO{counter,5} = mean(sum(counts(1,:,:,4),3));
%     significantGO{counter,6} = mean(sum(counts(1,:,:,5),3));
%     significantGO{counter,7} = std(sum(counts(1,:,:,1),3));
%     significantGO{counter,8} = std(sum(counts(1,:,:,2),3));
%     significantGO{counter,9} = std(sum(counts(1,:,:,3),3));
%     significantGO{counter,10} = std(sum(counts(1,:,:,4),3));
%     significantGO{counter,11} = std(sum(counts(1,:,:,5),3));
%     significantGO{counter,12} = pall_sorted(counter);
%     significantGO{counter,13} = sortedq(counter);
%     counter = counter + 1;
% end
% fileID = fopen('significantGOtitlesQvalsCorrectedGut.csv','w');
% formatSpec0 = '%s , %s , %s , %s , %s , %s, %s, %s, %s, %s, %s, %s, %s\n';
% header = {'GO Title','cecum count','ileum count','jejunum count','prox colon count','Stomach count','cecum std','ileum std','jejunum std','prox colon std','Stomach std','p-val','q-val'};
% fprintf(fileID,formatSpec0,header{1,:});
% formatSpec = '%s,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n';
% [nrows, ncols] = size(significantGO);
% for row = 1:nrows
%     fprintf(fileID,formatSpec,significantGO{row,:});
% end
% fclose(fileID);