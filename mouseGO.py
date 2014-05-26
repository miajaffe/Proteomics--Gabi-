# This Python script reads in the information from the goa file for mouse.
# The script then parses for only the proteins of interest (those in OverlordMatrix)
# and saves the GO codes and types (biological, molecular, cellular) to file.
import csv
# http://bmi214.stanford.edu/index.php?n=Site.PythonTutorial#toc18
# http://stackoverflow.com/questions/7856296/parsing-csv-tab-delimited-txt-file-with-python
allData = list(csv.reader(open('gene_association.goa_mouse','rb'),delimiter='\t'))
allData = allData[9:]
counter = 0
# GO code
value1 = []
# Biological, molecular, cellular
value2 = []
value = []
dataDic = {}
currentKey = allData[0][1]
for i in xrange(len(allData)):
    currentUniProt = allData[i][1]
    if currentUniProt == currentKey:
        if (len(value1) >= 1 and value1[len(value1)-1] != allData[i][4]):
            value1.append(allData[i][4])
            value2.append(allData[i][8])
    else:
        value = [value1,value2]
        dataDic[allData[i-1][1]] = value
        value = []
        value1 = []
        value2 = []
        currentKey = allData[i][1]
        value1.append(allData[i][4])
        value2.append(allData[i][8])
dataDic[currentUniProt] = value
axes = list(csv.reader(open('axes.txt','rb'),delimiter='\t'));
finalDic = {}
f = open('final_Dictionary_GO.txt','w')
for i in xrange(len(axes[0])):
    currentProtein = axes[0][i][0:7]
    print currentProtein
    if dataDic.has_key(currentProtein):
        # Adjust for MATLAB 1 indexing
        f.write(str(i+1) + '\n')
        f.write(str(dataDic[currentProtein][0]) + '\n')
        f.write(str(dataDic[currentProtein][1]) + '\n')
        f.write('ENDLINE' + '\n')
f.close()


    
    
