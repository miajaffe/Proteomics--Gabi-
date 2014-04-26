% This function returns the OverlordMatrix indeces associated with a given
% set of parameters.  These parameters can either be 1) the UniProt ID and
% the peptide letter sample name or 2) [UniProt ID, mouse ID, colonization state,
% location].  For example, ['Q8C5B4' 'Peptides A'] is equivalent to ['Q8C5B4'
% 'mouse3' 'BT' 'prox colon'].  The function takes either 4 or 6 inputs in
% the following possible orders: 1) [axes cell array, PeptideMap Map,
% UniProt ID, peptide letter sample name] or 2) [axes cell array, Peptide
% Map map, UniProt ID, mouse ID, colonization state,location along GI
% tract].  The output is always a 1 x 4 vector containing the indeces
% associated with the given conditions to be plugged into the
% OverlordMatrix to determine count value.
function outputIndeces = OverlordMatrixIndexGenerator(axes,PeptideMap,var1,var2,var3,var4)
% Axis 1 = UniProt ID
axis1key = axes{1};
% Axis 2 = mouse ID
axis2key = axes{2};
% Axis 3 = colonization state
axis3key = axes{3};
% Axis 4 = location along GI tract
axis4key = axes{4};
% Checks if either the UniProt ID and peptide letter sample name are
% given or the full set of [UniProt ID, mouse ID, colonization state,
% location] are given and checks the components of 'axes' to find the
% numerical value associated with each non-numerical value.
if nargin == 4
    index1 = find(ismember(axis1key,var1));
    peptide_data = PeptideMap(var2);
    index2 = peptide_data(1);
    index3 = peptide_data(2);
    index4 = peptide_data(3);
    outputIndeces = [index1,index2,index3,index4];
else
    index1 = find(ismember(axis1key,var1));
    index2 = find(ismember(axis2key,var2));
    index3 = find(ismember(axis3key,var3));
    index4 = find(ismember(axis4key,var4));
    outputIndeces = [index1,index2,index3,index4];
end
end