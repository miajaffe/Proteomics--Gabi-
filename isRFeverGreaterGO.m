% Checks to see if there are proteins where RF has greater counts than
% either BT or GF
load GOenrichMat
load IndextoGOConverterStr
load allGODic
RF_GF = GOenrichMat(:,:,3,:) - GOenrichMat(:,:,1,:);
RF_BT = GOenrichMat(:,:,2,:) - GOenrichMat(:,:,1,:);
RF_GF = mean(RF_GF,2);
RF_BT = mean(RF_BT,2);
RF_GF = sum(RF_GF,4);
RF_BT = sum(RF_BT,4);
for i = 1:1:2991
    if RF_GF(i) > 0
        temp = IndextoGOConverterStr(num2str(i));
        value = allGODic(temp);
        fprintf('GF')
        value{1}
    elseif RF_BT(i) > 0
        temp = IndextoGOConverterStr(num2str(i));
        value = allGODic(temp);
        fprintf('BT')
        value{1}
    end
end