 %Fecal data preparation


% Reads in excel file with data
%[fecalNum,fecalTxt]=xlsread('fecalMouseCoreProteome.xlsx');
% saved as mat files

% Initialization
% Load fecalNum.mat and fecalTxt.mat
fecalProtIds=fecalTxt(3:end,1);
fecalSpectralCounts=fecalNum(:,10);
% removes all proteins with spectral counts less than 5
% Adviced to do so by Josh as spectral counts below 5 are noise
notableProteinCount=find(fecalSpectralCounts > 5);
notableProteinIds=fecalProtIds(notableProteinCount);
% removed repeats of same protIds
coreFecalProteome=unique(notableProteinIds);
% removed decoys from dataset
coreFecalProteome=coreFecalProteome(31:end,1);

% shortens protein Id to enable comparison with gut protein IDs
% example shortens sp|A1L314|MPEG1_MOUSE to 'A1L314'
for i=1:length(coreFecalProteome)
    protId=coreFecalProteome{i};
    % takes substring of long protein id , 
    % example shortens sp|A1L314|MPEG1_MOUSE to 'A1L314'
    protId=protId(4:9);
    coreFecalProteome{i}=protId;
    
end