clear all
close all
clc
load('axes')
load('GOenrichMat')
load('GOtoIndexConverter')
load('IndextoGOConverter')
%% Let's see what differences there are between colonization states
GF_GO = sum(sum(GOenrichMat(:,:,1,:),2),4);
BT_GO = sum(sum(GOenrichMat(:,:,2,:),2),4);
RF_GO = sum(sum(GOenrichMat(:,:,3,:),2),4);
index = 1:1:length(GF_GO);
figure
scatter(index,GF_GO,'.r')
hold on
scatter(index,BT_GO,'.g');
scatter(index,RF_GO,'.b');
legend('GF','BT','RF')
xlabel('GO ID index')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index')
%% Let's see what differences there are between mice
maus1 = sum(sum(GOenrichMat(:,1,:,:),3),4);
maus2 = sum(sum(GOenrichMat(:,2,:,:),3),4);
maus3 = sum(sum(GOenrichMat(:,3,:,:),3),4);
index = 1:1:length(maus1);
figure
scatter(index,maus1,'.r')
hold on
scatter(index,maus2,'.g');
scatter(index,maus3,'.b');
legend('mouse1','mouse2','mouse3')
xlabel('GO ID index')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index')
%% Let's see what differences there are between locations
cecum = sum(sum(GOenrichMat(:,:,:,1),2),3);
ileum = sum(sum(GOenrichMat(:,:,:,2),2),3);
jejunum = sum(sum(GOenrichMat(:,:,:,3),2),3);
prox_colon = sum(sum(GOenrichMat(:,:,:,4),2),3);
Stomach = sum(sum(GOenrichMat(:,:,:,5),2),3);
figure
scatter(index,cecum,'.r')
hold on
scatter(index,ileum,'.g')
scatter(index,jejunum,'.b')
scatter(index,prox_colon,'.m')
scatter(index,Stomach,'.k')
legend('cecum','ileum','jejunum','prox colon','stomach')
xlabel('GO ID index')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index')
%% Let's calculate the standard deviation for each case at each GO ID and the mean of the std
%% Between colonization states
allCol = sum(sum(GOenrichMat(:,:,:,:),2),4);
allColStd = zeros(length(index),1);
normColStd = allColStd;
meanCol = allColStd;
for ii = 1:1:length(index)
    currStdVec = [allCol(ii,1,1,1),allCol(ii,1,2,1),allCol(ii,1,3,1)];
    allColStd(ii) = std(currStdVec);
    meanCol(ii) = mean(currStdVec);
    normColStd(ii) = std(currStdVec)/mean(currStdVec);
end
%% Between mice
allMaus = sum(sum(GOenrichMat(:,:,:,:),3),4);
allMausStd = zeros(length(index),1);
normMausStd = allMausStd;
meanMaus = allMausStd;
for ii = 1:1:length(index)
    currStdVec = [allMaus(ii,1,1,1),allMaus(ii,2,1,1),allMaus(ii,3,1,1)];
    allMausStd(ii) = std(currStdVec);
    meanMaus(ii) = mean(currStdVec);
    normMausStd(ii) = std(currStdVec)/mean(currStdVec);
end
%% Between locations
allLoc = sum(sum(GOenrichMat(:,:,:,:),2),3);
allLocStd = zeros(length(index),1);
normLocStd = allLocStd;
meanLoc = allLocStd;
for ii = 1:1:length(index)
    currStdVec = [allLoc(ii,1,1,1),allLoc(ii,1,1,2),allLoc(ii,1,1,3),allLoc(ii,1,1,4),allLoc(ii,1,1,5)];
    allLocStd(ii) = std(currStdVec);
    normLocStd(ii) = std(currStdVec)/mean(currStdVec);
    meanLoc(ii) = mean(currStdVec);
end
%%
avgStdCol = mean(allColStd);
avgStdMaus = mean(allMausStd);
avgStdLoc = mean(allLocStd);
fprintf('The avg std for colonization state is %f\n',avgStdCol)
fprintf('The avg std for mouse replicated is %f\n',avgStdMaus)
fprintf('The avg std for location is %f\n',avgStdLoc)
%% An interesting question: which GO IDs have greatest % variance?
%% Between colonization states
% Look for the 10 most frequent GO codes with mean normalized counts >= 0.5
[sortedColStdValue,sortedColStd] = sort(normColStd,'descend');
meanCol = meanCol(sortedColStd);
count = 0;
count2 = 1;
threshold = 0.5;
cutoff = 20;
highFreqCol = {};
xaxis = 1:1:cutoff;
highFreqColMat = zeros(cutoff,3);
while count < cutoff || count2 > length(meanCol)
    if meanCol(count2) >= threshold
        count = count + 1;
        tempIndex = num2str(sortedColStd(count2));
        highFreqCol{count} = IndextoGOConverter(tempIndex);
        highFreqColMat(count,1) = GF_GO(sortedColStd(count2));
        highFreqColMat(count,2) = BT_GO(sortedColStd(count2));
        highFreqColMat(count,3) = RF_GO(sortedColStd(count2));
    end
    count2 = count2 + 1;
end
figure
scatter(xaxis,highFreqColMat(:,1),'r')
hold on
scatter(xaxis,highFreqColMat(:,2),'g')
scatter(xaxis,highFreqColMat(:,3),'b')
legend('GF','BT','RF')
xlabel('GO ID (not index)')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index: Greatest Reasonable % Variation')
set(gca,'fontsize',8,'XLim',[1 cutoff],'XTick',[1:1:cutoff],'XTickLabel',highFreqCol)
%% Between Mice
[sortedMausStdValue,sortedMausStd] = sort(normMausStd,'descend');
meanMaus = meanMaus(sortedMausStd);
count = 0;
count2 = 1;
threshold = 0.5;
cutoff = 20;
highFreqMaus = {};
xaxis = 1:1:cutoff;
highFreqMausMat = zeros(cutoff,3);
while count < cutoff || count2 > length(meanMaus)
    if meanMaus(count2) >= threshold
        count = count + 1;
        tempIndex = num2str(sortedMausStd(count2));
        highFreqMaus{count} = IndextoGOConverter(tempIndex);
        highFreqMausMat(count,1) = maus1(sortedMausStd(count2));
        highFreqMausMat(count,2) = maus2(sortedMausStd(count2));
        highFreqMausMat(count,3) = maus3(sortedMausStd(count2));
    end
    count2 = count2 + 1;
end
figure
scatter(xaxis,highFreqMausMat(:,1),'r')
hold on
scatter(xaxis,highFreqMausMat(:,2),'g')
scatter(xaxis,highFreqMausMat(:,3),'b')
legend('Mouse1','Mouse2','Mouse3')
xlabel('GO ID (not index)')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index: Greatest Reasonable % Variation')
set(gca,'XTickLabel',highFreqMaus,'fontsize',8,'XLim',[1 cutoff],'XTick',[1:1:cutoff])

%% Between locations
[sortedLocStdValue,sortedLocStd] = sort(normLocStd,'descend');
meanLoc = meanLoc(sortedLocStd);
count = 0;
count2 = 1;
threshold = 0.5;
cutoff = 20;
highFreqLoc = {};
xaxis = 1:1:cutoff;
highFreqLocMat = zeros(cutoff,3);
while count < cutoff || count2 > length(meanLoc)
    if meanLoc(count2) >= threshold
        count = count + 1;
        tempIndex = num2str(sortedLocStd(count2));
        highFreqLoc{count} = IndextoGOConverter(tempIndex);
        highFreqLocMat(count,1) = cecum(sortedLocStd(count2));
        highFreqLocMat(count,2) = ileum(sortedLocStd(count2));
        highFreqLocMat(count,3) = jejunum(sortedLocStd(count2));
        highFreqLocMat(count,4) = prox_colon(sortedLocStd(count2));
        highFreqLocMat(count,5) = Stomach(sortedLocStd(count2));
    end
    count2 = count2 + 1;
end
figure
scatter(xaxis,highFreqLocMat(:,1),'r')
hold on
scatter(xaxis,highFreqLocMat(:,2),'g')
scatter(xaxis,highFreqLocMat(:,3),'b')
scatter(xaxis,highFreqLocMat(:,4),'m')
scatter(xaxis,highFreqLocMat(:,5),'k')
legend('cecum','ileum','jejunum','prox colon','stomach')
xlabel('GO ID (not index)')
ylabel('Normalized Counts')
title('Normalized Counts vs. GO ID index: Greatest Reasonable % Variation')
set(gca,'XTickLabel',highFreqLoc,'fontsize',8,'XLim',[1 cutoff],'XTick',[1:1:cutoff])
%% Compile high variants
allVariants = [highFreqCol;highFreqMaus;highFreqLoc];
allDefinitions = {};
counter = 0;
for ii = 1:1:size(allVariants,1)
    for iii = 1:1:size(allVariants,2)
        counter = counter + 1
        tempGO = num2str(allVariants{ii,iii});
        lenZeros = 7 - length(tempGO);
        tempZeros = '';
        if lenZeros >0
            for j=1:1:lenZeros
                tempZeros = strcat(tempZeros,'0');
            end
        end
        allVariants{ii,iii} = strcat('GO:',tempZeros,tempGO);
        urlStr = allVariants{ii,iii};
        url = strcat('http://amigo.geneontology.org/amigo/term/',urlStr);
        tempData = urlread(url);
        startI = findstr(tempData,'>Definition');
        endI = findstr(tempData,'<em>');
        roughData = tempData(startI:endI(2));
        polishedData = roughData(49:length(roughData)-5);
        allDefinitions{ii,iii} = polishedData;
    end
end
%%
fileID1 = fopen('InterestingGOIDs.dat','w');
fileID2 = fopen('InterestingDefinitions.dat','w');
formatSpec = '%s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n %s\n\n\n\n\n\n\n\n\n\n\n';
[nrows,ncols] = size(allVariants);
for row = 1:nrows
    fprintf(fileID1,formatSpec,allVariants{row,:});
    fprintf(fileID2,formatSpec,allDefinitions{row,:});
end
fclose(fileID1);
fclose(fileID2);

%% Find definitions of each GO Term of interest

% set(gca,'XTickLabel',highFreqLoc,'fontsize',8)
% %% PCA <-- Fuck this, getting lost in the maths, not enough intuition
% dims = [size(GOenrichMat,1),size(GOenrichMat,2),size(GOenrichMat,3),size(GOenrichMat,4)];
% % Prep matrix for covariance matrix is 5D: 
% % (Value,GO ID,Mouse,Colonization,Location)
% % prepMat = zeros(dims(1),dims(2)*dims(3)*dims(4));
% prepMat = [];
% counter = 0;
% % Taken from PCA_andHierarchicalClustering.m
% all_labels = {};
%     for iii = 1:1:dims(2)
%         for iv = 1:1:dims(3)
%             for ivi = 1:1:dims(4)
%                 prepMat = [prepMat GOenrichMat(:,iii, iv, ivi)];
%                 label = strcat(axes{2}{iii}, '_', axes{3}{iv} , '_', axes{4}{ivi});           
%                 all_labels = [all_labels label];
%             end
%         end
%     end
% [pc,score,latent,tsquare] = princomp(prepMat);
% cumsum(latent)./sum(latent);
% figure
% biplot(pc(:,1:2),'Scores',score(:,1:2),'VarLabels',all_labels)
% % allAvg = mean(prepMat);
% % nv = size(prepMat,1);
% % diff_avg = prepMat - repmat(allAvg,nv,1);
% % ndim = size(prepMat,2);
% % cov = zeros(ndim,ndim);
% % for ii = 1:ndim
% %     for iii = 1:ndim
% %         cov(ii,iii) = 1/(nv-1)*sum(diff_avg(:,ii).*diff_avg(:,iii));
% %     end
% % end
% % figure
% % imagesc(cov);
% % colorbar
% % [V,D] = eigs(cov);
% % eigenvalues = diag(D);
% % variances = eigenvalues/sum(eigenvalues);
% % for i=1:5
% %  fprintf('The fraction of the variance from mode %d is %g\n',i,variances(i));
% % end
% % %%
% % % step 4: calculate eigs of covariance_matrix
% % 
% % [eigenvectors,eigenvalues,explained] = pcacov(cov);
% % 
% % 
% % %% 
% % % step 5: generate files for a PCA plot
% % 
% % eigenvector_matrix = eigenvectors;
% % add_labels = {'eigvals'; '% variation explained'};
% % add_title = {'pc vector number'};
% % pc_vector_labels = all_labels';
% % pc_vector_number = vertcat(add_title,pc_vector_labels, add_labels);
% % 
% % pc_labels = (1:size(all_labels,2));
% % eigvals = eigenvalues';
% % variation_explained = explained';
% % 
% % eigenvector_matrix_labeled = vertcat(pc_labels, eigenvector_matrix, eigvals, variation_explained);
% % eigenvector_table = table(pc_vector_number, eigenvector_matrix_labeled);
% % 
% % writetable(eigenvector_table,'QIIME_proteomics_pc.txt','Delimiter','\t');