from urllib import request

from porestat.utils.DataFrame import DataFrame

allgenomes = DataFrame.parseFromFile("../../../ena_bacteria_list.csv")
print(allgenomes.getHeader())

for row in allgenomes:

    protInfo = row['proteins']
    print(protInfo)


for row in allgenomes:

    protInfo = row['proteins']

    if protInfo==None or len(protInfo) == 0 or len(protInfo.strip()) == 0 or protInfo == 'n/a' or protInfo == 'None':
        continue

    downloadFile = row['seqID'] + ".gb"
    downloadLocation = "../../../genomes/"
    print(downloadFile)

    request.urlretrieve("http://www.ebi.ac.uk/ena/data/view/"+row['seqID']+"&display=txt&expanded=true", downloadLocation + "/" + downloadFile)





