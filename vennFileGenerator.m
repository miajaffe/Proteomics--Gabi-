%%Prints out tab delimited txt file of the proteins present in the different 
%locations in the GI tract

% Initialization
load axes.mat, load normOverlord_shannon

protIds=axes{1};

% initialize 
% 1 = cecum 2=ileum 3=jejunum 4=proxCol 5=stomach

%find unique protIds from all replicates in each location
indicesOF1ColonizationState=cell(5,3);

uniqueInd=cell(1,5); % unique Indices of proteins of all replicates in a particular colonization 
%state and GI location. For instance the unique indices of all replicates
%in the cecum in Germ free colonization state
uniqueProtIds=cell(3,5); % unique protIds of all replicates in every colonization state
% and GI location
for GIlocation=1:5
    for colonisationState=1:3
        for mouseReplicate=1:3
            protAbundanceInReplicate=normOverlord_shannon(:,mouseReplicate,colonisationState,GIlocation);
            indicesOF1ColonizationState{GIlocation,mouseReplicate}=find(protAbundanceInReplicate);   
        end
        %saves unique indices of all mouse replicates in a given
        %colonization state in a cell that will contain all unique
        concatInd=cat(1,indicesOF1ColonizationState{GIlocation,1},indicesOF1ColonizationState{GIlocation,2},indicesOF1ColonizationState{GIlocation,3});
        uniqueInd{1,    GIlocation}=unique(concatInd);
        %finds associated protIds of the unique indices in a given location
        Id=cell(1,length(uniqueInd{1,GIlocation}));
        indices=uniqueInd{1,GIlocation}; % unique indices of protein Ids
        % Matching unique indices with prot Ids
        for i=1:length(indices)
   
            Id{1,i}=protIds{1,indices(i,1)};
            
        end
        uniqueProtIds{colonisationState,GIlocation}=Id;
    end
    
    
    
end

