% Mia GO-code processing
clear all
close all
clc
[num,txt,raw] = xlsread('protein_IDs_from_clusters.xlsx','A1:R42');
load axes
load GOArray

%% B-theta
counter = 3;
for i = 1:1:5
    while length(txt{counter,i}) > 0
        BTGO{counter-2,i} = txt{counter,i};
        counter = counter + 1;
    end
    counter = 3;
end
%% Germ-Free
counter = 3;
for i = 7:1:11
    while counter < 42 && length(txt{counter+1,i}) > 0
        GFGO{counter-2,i-6} = txt{counter,i};
        counter = counter + 1;
    end
    if counter == 42
        GFGO{counter-2,i-6} = txt{counter,i};
    end
    counter = 3;
end
%% Conventionally raised
counter = 3;
for i = 13:1:18
    while length(txt{counter,i}) > 0
        RFGO{counter-2,i-12} = txt{counter,i};
        counter = counter + 1;
    end
    counter = 3;
end

%% Aggregate GO terms for each case
% For each colonization state, we make a cell array of Maps.  The idea is
% that each cell in the cell array corresponds to one cluster in a given
% colonization state, which in turn corresponds to a given map.  The kys
% give a non-redundant list of all occuring GO codes and the values give a
% tally for the abundance for each GO code.
%% B-theta
counter = 1;
BTGOMaps = {};
for i = 1:1:size(BTGO,2)
    BTGOMaps{i} = containers.Map();
    while counter <= size(BTGO,1) && length(BTGO{counter,i}) > 0
        protindex = find(ismember(axes{1},BTGO{counter,i}));
        currGOIDs = GOArray{protindex,1};
        for j = 1:1:length(currGOIDs)
            if isKey(BTGOMaps{i},currGOIDs{j})
                value = BTGOMaps{i}(currGOIDs{j}) + 1;
                BTGOMaps{i}(currGOIDs{j}) = value;
            else
                BTGOMaps{i}(currGOIDs{j}) = 1;
            end
        end
        counter = counter + 1;
    end
    counter = 1;
end
%% GF
counter = 1;
GFGOMaps = {};
for i = 1:1:size(GFGO,2)
    GFGOMaps{i} = containers.Map();
    while counter <= size(GFGO,1) && length(GFGO{counter,i}) > 0
        protindex = find(ismember(axes{1},GFGO{counter,i}));
        currGOIDs = GOArray{protindex,1};
        for j = 1:1:length(currGOIDs)
            if isKey(GFGOMaps{i},currGOIDs{j})
                value = GFGOMaps{i}(currGOIDs{j}) + 1;
                GFGOMaps{i}(currGOIDs{j}) = value;
            else
                GFGOMaps{i}(currGOIDs{j}) = 1;
            end
        end
        counter = counter + 1;
    end
    counter = 1;
end
%% RF
counter = 1;
RFGOMaps = {};
for i = 1:1:size(RFGO,2)
    RFGOMaps{i} = containers.Map();
    while counter <= size(RFGO,1) && length(RFGO{counter,i}) > 0
        protindex = find(ismember(axes{1},RFGO{counter,i}));
        currGOIDs = GOArray{protindex,1};
        for j = 1:1:length(currGOIDs)
            if isKey(RFGOMaps{i},currGOIDs{j})
                value = RFGOMaps{i}(currGOIDs{j}) + 1;
                RFGOMaps{i}(currGOIDs{j}) = value;
            else
                RFGOMaps{i}(currGOIDs{j}) = 1;
            end
        end
        counter = counter + 1;
    end
    counter = 1;
end
%% Analysis of keys, i.e. unique GO codes
% http://stackoverflow.com/questions/12415361/intersection-of-2-cell-arrays-of-strings-in-matlab
% Let's compare GF{4} and BT{4} since they have strongest K-means
% association:
GF4keys = keys(GFGOMaps{4});
BT4keys = keys(BTGOMaps{4});
GF4val = values(GFGOMaps{4});
GF4sum = 0;
for i = 1:1:length(GF4val)
    GF4sum = GF4sum + GF4val{i};
end
BT4val = values(BTGOMaps{4});
BT4sum = 0;
for i = 1:1:length(BT4val)
    BT4sum = BT4sum + BT4val{i};
end
intrsctGOIDs = GF4keys(ismember(GF4keys,BT4keys));
intrsctsumGF = 0;
intrsctsumBT = 0;
for i = 1:1:length(intrsctGOIDs)
    intrsctsumGF = intrsctsumGF + GFGOMaps{4}(intrsctGOIDs{i});
    intrsctsumBT = intrsctsumBT + BTGOMaps{4}(intrsctGOIDs{i});
end
GFratio = intrsctsumGF/GF4sum
BTratio = intrsctsumBT/BT4sum