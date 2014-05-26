load('axes')
fileID = fopen('axes.txt','w');
formatSpec = '%s\t';
for i = 1:1:size(axes{1},2)
    fprintf(fileID,formatSpec,axes{1}{i});
end
fclose(fileID);