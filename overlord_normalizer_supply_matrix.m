% Outputs the normalized OverlordMatrix normalized both across rows and
% across columns.  Normalization by columns (down every row) allows
% comparisons between relative abundance of proteins in a given sample,
% e.g. we can compare protein cadherin-1 to protein IL-2 within sample
% Peptide A.
% Normalization by rows (across every column), i.e. comparisons can be made
% for the relative abundance of a given protein across all samples, e.g. we
% can "compare" protein cadherin-1's relative abundance across all samples
% (Peptides A to ZZ).  Also required is the MOUSE.fasta downloaded from
% UniProt!  Link:
% http://www.uniprot.org/uniprot/?query=taxonomy%3a10090&sort=score&format=*
% hit download
% UPDATE: If you have MusProt.mat in the working directory, the FASTA data
% is not needed and the program runs fastA.


% [redunancies] = the vector of indeces within axes for redundant protein
% ID's.

% Form to use if MusProtRaw.mat unavailable
% function [normOverlord] = OverlordNormalizer(choice,MusProtRaw)
% Form to use if MusProtRaw.mat is available
function [normOverlord,redundancies] = overlord_normalizer_supply_matrix()
% Normalized by columns, compare all proteins within one sample.  To do
% this, we must divide each protein's number of counts by the protein amino
% acid length, since number of counts scales with number of enzymatic
% cleavage sites, which in turn scales with the length of the protein.  For
% future versions (not in the context of this class), a more precise
% normalization is to divide by the number of specific cleavage sites used
% for the sample processing (e.g. only number of trypsin sites if trypsin
% is used)
% If MusProtRaw.mat is unavailable
% MusProt = fastaread(MusProtRaw);
% If MusProtRaw.mat is available
load('MusProt.mat')
load('OverlordMatrix.mat');
load('axes.mat');

%%
cutoff = 5;
cutoff_matrix = OverlordMatrix;
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            [rows] = find(OverlordMatrix(:,mouse_num, colonization, loc) <= cutoff);
            cutoff_matrix(rows,mouse_num, colonization, loc) = 0;
        end
    end
end
OverlordMatrix2 = cutoff_matrix
%%
% First step is to acquire the MW of each and every protein.
proteinID = {MusProt.Header};
proteinSeq = {MusProt.Sequence};
% Key is protein ID, value is length of protein in # amino acids
lengthMap = containers.Map();
for ii = 1:1:size(proteinID,2)
    temp = proteinID{ii};
    temp = temp(4:9);
    lengthSeq = length(proteinSeq{1,ii});
    lengthMap(temp) = lengthSeq;
end
colNormFactor = zeros([1 size(axes{1,1},2)]);
keysProt = axes{1,1};
noGood = {};
ggCounter = 0;
redundancies = [];
for ii = 1:1:length(colNormFactor)
    %         fprintf('%d\n',ii)
    % Ignores decoys with '_' in ID name
    if length(strfind(keysProt{ii}, '_')) == 0
        if isKey(lengthMap,keysProt{ii}) == 1
            colNormFactor(ii) = lengthMap(keysProt{ii});
        else
            % To deal with weird reduncancies, search the web to find
            % AA length of equivalent protein
            % http://stackoverflow.com/questions/8061344/how-to-search-for-a-string-in-cell-array-in-matlab
            redundancies(length(redundancies)+1) = find(ismember(axes{1},keysProt{ii}));
            url = strcat('http://www.uniprot.org/uniprot/',keysProt{ii});
            urldata = urlread(url);
            urlposition = findstr(urldata,' AA');
            strcount = urlposition-1;
            while isstrprop(urldata(strcount), 'digit') == 1
                strcount = strcount - 1;
            end
            templength = str2num(urldata(strcount+1:urlposition));
            if length(urlposition) > 0
                colNormFactor(ii) = templength;
                % To deal with deleted UniProt entries, e.g.D3YUL9
            else
                colNormFactor(ii) = 1;
            end
        end
    else
        % Prevent div by 0 downstream
        colNormFactor(ii) = 1;
    end
end
% Tile the colNormFacx
tiledNorm = repmat(colNormFactor',1,size(OverlordMatrix,2),size(OverlordMatrix,3),size(OverlordMatrix,4));
normOverlord = OverlordMatrix ./ tiledNorm;
% Normalized by rows, compare a given protein within all samples.  To do
% this, we must divide every protein's number of counts by the total number
% of counts within= the given sample.  This is because, even though the
% instrument and methods are generally the same, small differences affect
% the overall responsiveness to a given sample.


% Normalize across samples
% Find the sum of all counts in every sample
rowNormFactor = sum(normOverlord,1);
% Tile the rowNormFactor so that it has the same dimensions as
% normOverlord
tiledNorm = repmat(rowNormFactor,size(normOverlord,1),1,1,1);
% Divide normOverlord by the tiledNorm
normOverlord_shannon = normOverlord ./ tiledNorm;
save('normOverlord_shannon','normOverlord_shannon')

end