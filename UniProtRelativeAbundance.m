% Examines incidence of proteins of interest within
% normOverlordFinal 
clear all
close all
load axes
load GOtoIndexConverterStr
load IndextoGOConverterStr
load allGODic
load titleToGO
load normOverlordFinal
%% Find GO codes with q-value <= arbitrary such that less than 1 significant sample is a false positive
allProt = axes{1};
for i = 1:1:length(allProt)
    Protcurr = allProt{i};
    for j = 1:1:3
        GF_BT_RF_Prot1(j,1) = sum(normOverlordFinal(i,j,1,:),4);
        GF_BT_RF_Prot1(j,2) = sum(normOverlordFinal(i,j,2,:),4);
        GF_BT_RF_Prot1(j,3) = sum(normOverlordFinal(i,j,3,:),4);
    end
    [p,table,stats] = anova1(GF_BT_RF_Prot1,{'GF','BT','RF'},'off');
    pall(i) = p;
    statsall(i) = stats;
end
pall_sorted = sort(pall,'ascend');
[FDR, q] = mafdr(pall);
[sortedq, qind] = sort(q,'ascend');
counter = 1;
qcutoff = 0.03;
statssorted = statsall(qind);
while sortedq(counter) <= qcutoff
    index = qind(counter);
    temp = axes{1}{index};
    counts = normOverlordFinal(index,:,:,:);
    significantGO{counter,1} = temp;
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
fileID = fopen('significantUniProtQvalsCorrected5pctFDR.csv','w');
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
