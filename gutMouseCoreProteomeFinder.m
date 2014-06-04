% Finds the core proteome in each location across all colonisation states

% Initialization
% load axes.mat and normOverlord_shannon.mat

protIds=axes{1};

% initialize
% 1 = cecum 2=ileum 3=jejunum 4=proxCol 5=stomach

%find unique protIds from all replicates in each location
indicesOF1ColonizationState=cell(5,3);
indicesOFAllColonizationState=cell(1,3);

uniqueInd=cell(1,5); % unique Indices of proteins of all replicates in a particular colonization
%state and GI location. For instance the unique indices of all replicates
%in the cecum in Germ free colonization state
uniqueProtIds=cell(3,5); % unique protIds of all replicates in every colonization state
% and GI location

for colonisationState=1:3
    
    for GIlocation=1:5
        
        %         disp('this is the col')
        %         disp(colonisationState)
        for mouseReplicate=1:3
            protAbundanceInReplicate=normOverlord_shannon(:,mouseReplicate,colonisationState,GIlocation);
            
            indicesOF1ColonizationState{GIlocation,mouseReplicate}=find(protAbundanceInReplicate);
            
            
        end
        
    end
    indicesOFAllColonizationState{1,colonisationState}=indicesOF1ColonizationState;
    
    
end
germFree=indicesOFAllColonizationState{1};
bTheta=indicesOFAllColonizationState{2};
conventional=indicesOFAllColonizationState{3};
allProtIndicesPerColonisationState=cell(5,3);
for i=1:5
    
        % put all protIds in one cell
       % allProtIdsPerColonisationState{1,i}=cat(1,germFree{i,1},germFree{i,2},germFree{i,3},bTheta{i,1},bTheta{i,2},bTheta{i,3},conventional{i,1},conventional{i,2},conventional{i,3});
       allProtIndicesPerColonisationState{i,1}=cat(1,germFree{i,1},germFree{i,2},germFree{i,3});
       allProtIndicesPerColonisationState{i,2}=cat(1,bTheta{i,1},bTheta{i,2},bTheta{i,3});
       allProtIndicesPerColonisationState{i,3}=cat(1,conventional{i,1},conventional{i,2},conventional{i,3});
    
end 

% Find the intersection of all proteins indices in each colonization state
% to find core proteome per location. Intersection finds proteins that are
% found in each colonization state. This proteins are our core proteome of
% that location. 
intersectionOfProtIndicesPerColonisationState=cell(1,5);

%Find unique proteins in each colonization state
% uniqueProtIdsPerColonisationState=cell(1,5);
% 
for i=1:5
    gFBtIntersect=intersect(allProtIndicesPerColonisationState{i,1},allProtIndicesPerColonisationState{i,2});
    allIntersect=intersect(gFBtIntersect,allProtIndicesPerColonisationState{i,3});
    intersectionOfProtIndicesPerColonisationState{1,i}=allIntersect;
    %intersectionOfProtIdsPerColonisationState{1,i}=unique(allProtIdsPerColonisationState{1,i});
end 

coreGutProteomePerGILocation=cell(1,5);
for i=1:length(intersectionOfProtIndicesPerColonisationState)
    % fill intersectionProteinIdsPercolonisationstate with all protIds
    % associated with indices in the intersected proteins. This protein ids
    % form the core proteome at the given GI location
    coreGutProteomePerGILocation{1,i}=protIds(intersectionOfProtIndicesPerColonisationState{i});
end
