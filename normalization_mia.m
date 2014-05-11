% This function normalizes the values in the Overlord matrix such that each
% count value is divided by the total protein counts in that particular
% sample, then mutliplied by a factor.
% The normlized matrix is saved in a new variable 'NewOverlord'

function NewOverlord = normalization_mia(OverlordMatrix)

% Loop through all three mouse replicates:
for i = 1:3
    % Loop through all three colonization states:
    for j = 1:3
       % Loop through all 5 locations:
        for k = 1:5 
           % Normalize protein counts by dividing by total of all counts and multiplying by a constant 
            NewOverlord(:,i,j,k) = OverlordMatrix(:,i, j, k)/sum(OverlordMatrix(:,i, j, k))*1000;
       end
    end
end

