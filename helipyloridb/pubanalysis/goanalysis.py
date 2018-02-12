from collections import defaultdict
from collections import Counter
import operator
from nertoolkit.geneontology.GeneOntology import GeneOntology

locusTag2GO = {'HP_1213': {'GO:0006396', 'GO:0003723', 'GO:0006402', 'GO:0004654', 'GO:0046872', 'GO:0005737'}, 'HP_1141': {'GO:0004479'}, 'HP_0091': {'GO:0009307', 'GO:0009036', 'GO:0003677'}, 'HP_0818': {'GO:0005886', 'GO:0005215', 'GO:0006810', 'GO:0016021'}, 'HP_1486': {'GO:0016021'}, 'HP_0770': {'GO:0009306', 'GO:0005886', 'GO:0044780', 'GO:0016021'}, 'HP_0415': {'GO:0005886', 'GO:0055085', 'GO:0016021'}, 'HP_1349': {'GO:0016021'}, 'HP_0860': {'GO:0034200', 'GO:0009244', 'GO:0000287', 'GO:0097171', 'GO:0008270', 'GO:0005737'}, 'HP_0048': {'GO:0008270', 'GO:0046944', 'GO:0003725', 'GO:0016787', 'GO:0016743'}, 'HP_0142': {'GO:0006284', 'GO:0003824', 'GO:0003677'}, 'HP_0312': {'GO:0005524'}, 'HP_1282': {'GO:0046872', 'GO:0000162', 'GO:0004049'}, 'HP_0272': {'GO:0016021'}, 'HP_0043': {'GO:0000271', 'GO:0016853', 'GO:0016779'}, 'HP_0731': {'GO:0005525'}, 'HP_0289': {'GO:0042802'}, 'HP_0661': {'GO:0004523', 'GO:0003676', 'GO:0005737', 'GO:0046872'}, 'HP_1499': {'GO:0003824'}, 'HP_0558': {'GO:0033817', 'GO:0006633'}, 'HP_0807': {'GO:0009279', 'GO:0004872', 'GO:0006810'}, 'HP_0793': {'GO:0042586', 'GO:0046872', 'GO:0006412'}, 'HP_1184': {'GO:0015297', 'GO:0015238', 'GO:0016021'}, 'HP_0462': {'GO:0006304', 'GO:0003677'}, 'HP_0547': {'GO:0019534'}, 'HP_0759': {'GO:0015297', 'GO:0015238', 'GO:0016021'}, 'HP_1072': {'GO:0004008', 'GO:0005886', 'GO:0005524', 'GO:0046872', 'GO:0016021'}, 'HP_0338': {'GO:0016021'}, 'HP_0656': {'GO:0051539', 'GO:0016765', 'GO:0046992', 'GO:0005506', 'GO:0009234'}, 'HP_0580': {'GO:0016021'}, 'HP_1100': {'GO:0051539', 'GO:0009255', 'GO:0046872', 'GO:0004456', 'GO:0019521'}, 'HP_0854': {'GO:0006163', 'GO:1902560', 'GO:0003920'}, 'HP_1471': {'GO:0006304', 'GO:0003677'}, 'HP_0342': {'GO:0016021'}, 'HP_1229': {'GO:0009088', 'GO:0004072', 'GO:0005524', 'GO:0019877', 'GO:0009089'}, 'HP_0599': {'GO:0016020', 'GO:0004871'}, 'HP_1566': {'GO:0016021'}, 'HP_1105': {'GO:0016757'}, 'HP_1044': {'GO:0046872', 'GO:0016787'}, 'HP_0608': {'GO:0016021'}, 'HP_0927': {'GO:0004222', 'GO:0005886', 'GO:0046872', 'GO:0016021'}, 'HP_0099': {'GO:0004871', 'GO:0016021'}, 'HP_0747': {'GO:0008176'}, 'HP_0743': {'GO:0007049', 'GO:0005886', 'GO:0071555', 'GO:0009252', 'GO:0008955', 'GO:0008360', 'GO:0051301', 'GO:0016021'}, 'HP_0995': {'GO:0006310', 'GO:0003677', 'GO:0015074'}, 'HP_0566': {'GO:0009089', 'GO:0005737', 'GO:0008837'}, 'HP_1361': {'GO:0016021'}, 'HP_0898': {'GO:0046872'}, 'HP_1259': {'GO:0070403', 'GO:0005737', 'GO:0036054', 'GO:0036055'}, 'HP_1046': {'GO:0042274', 'GO:0005737'}, 'HP_0755': {'GO:0008641'}, 'HP_0286': {'GO:0004222', 'GO:0005524', 'GO:0051301', 'GO:0016021'}, 'HP_0104': {'GO:0009166', 'GO:0046872', 'GO:0016788', 'GO:0000166'}, 'HP_1517': {'GO:0008168', 'GO:0006304', 'GO:0003677'}, 'HP_1277': {'GO:0004834'}, 'HP_0681': {'GO:0016021'}, 'HP_0235': {'GO:0005576', 'GO:0046677', 'GO:0008800'}, 'HP_1513': {'GO:0097056', 'GO:0001514', 'GO:0005737', 'GO:0004125'}, 'HP_0012': {'GO:0003896', 'GO:0008270', 'GO:0003677', 'GO:1990077'}, 'HP_1533': {'GO:0006235', 'GO:0006231', 'GO:0050797', 'GO:0050660'}, 'HP_0298': {'GO:0043190', 'GO:0055085'}, 'HP_0517': {'GO:0005886', 'GO:0005525', 'GO:0042254', 'GO:0005737', 'GO:0019843'}, 'HP_0238': {'GO:0002161', 'GO:0005524', 'GO:0006433', 'GO:0005737', 'GO:0004827'}}


count2locusTag = {}

count2locusTag[1] = ['HP_0012', 'HP_0030', 'HP_0036', 'HP_0043', 'HP_0047', 'HP_0059', 'HP_0064', 'HP_0091', 'HP_0099', 'HP_0104', 'HP_0142', 'HP_0235', 'HP_0238', 'HP_0262', 'HP_0286', 'HP_0289', 'HP_0298', 'HP_0312', 'HP_0338', 'HP_0342', 'HP_0415', 'HP_0421', 'HP_0428', 'HP_0456', 'HP_0465', 'HP_0517', 'HP_0519', 'HP_0547', 'HP_0558', 'HP_0566', 'HP_0568', 'HP_0580', 'HP_0599', 'HP_0608', 'HP_0681', 'HP_0743', 'HP_0747', 'HP_0755', 'HP_0759', 'HP_0762', 'HP_0770', 'HP_0783', 'HP_0788', 'HP_0793', 'HP_0807', 'HP_0818', 'HP_0854', 'HP_0860', 'HP_0868', 'HP_0884', 'HP_0898', 'HP_0906', 'HP_0927', 'HP_0995', 'HP_1044', 'HP_1046', 'HP_1072', 'HP_1078', 'HP_1100', 'HP_1107', 'HP_1141', 'HP_1157', 'HP_1177', 'HP_1184', 'HP_1213', 'HP_1229', 'HP_1259', 'HP_1265', 'HP_1277', 'HP_1282', 'HP_1361', 'HP_1471', 'HP_1486', 'HP_1499', 'HP_1513', 'HP_1517', 'HP_1566']
count2locusTag[2] = ['HP_0048', 'HP_0108', 'HP_0272', 'HP_0430', 'HP_0462', 'HP_0656', 'HP_0731', 'HP_0897', 'HP_1349', 'HP_1533']
count2locusTag[3] = ['HP_1105', 'HP_1455']
count2locusTag[5] = ['HP_0661']


goID2label = {}

allBaseGOs = set()
goCat2Count = Counter()

for x in locusTag2GO:
    for y in locusTag2GO[x]:
        allBaseGOs.add(y)
        goCat2Count[y] += 1

goObo = GeneOntology("/mnt/c/ownCloud/data/" + "miRExplore/go/go.obo")

locusTag2rootGO = defaultdict(set)

for locusTag in locusTag2GO:

    rootElems = set()

    testElems = list(locusTag2GO[locusTag])
    ignoreElems = set()

    allparents = set()

    bestElem = testElems[0]
    for elem in testElems:
        if goCat2Count[elem] > goCat2Count[bestElem]:
            bestElem = elem

    for elem in testElems:
        if goCat2Count[elem] == goCat2Count[bestElem]:
            rootElems.add(elem)

    testElems = []


    for i in range(0, len(testElems)):

        testElem = testElems[i]

        if testElem in rootElems:
            continue # should not happen

        if testElem in ignoreElems:
            continue

        res = goObo.getID(testElem)

        goID2label[res.id] = res.name

        if res == None:
            print("Not Found: ", res, testElem)

        setAllChildren = res.getAllChildren()

        if len(allparents) > 0:
            allparents = allparents.intersection(res.getAllParents())
        else:
            allparents = res.getAllParents()


        #print(len([x for x in setAllChildren if x.termid in testElems]))

        for child in setAllChildren:
            if child.termid in testElems[i:len(testElems)]:
                ignoreElems.add(child.id)

            if child.termid in rootElems:
                rootElems.remove(child.termid)

            if child.termid in testElems:
                print(testElems)

        parentInBase = False
        for parent in res.is_a:
            goID2label[parent.term.id] = parent.term.name

            if parent.termid in allBaseGOs:
                parentInBase = True


        if parentInBase:
            for parent in res.is_a:
                if parent.termid in allBaseGOs:
                    rootElems.add(parent.termid)

        else:
            rootElems.add(res.id)


    for goID in rootElems:
        res = goObo.getID(goID)
        goID2label[res.id] = res.name

    locusTag2rootGO[locusTag] = rootElems


import numpy as np
import matplotlib.pyplot as plt


for diff in count2locusTag:

    classCounts = Counter()

    for locusTag in count2locusTag[diff]:

        elems = locusTag2rootGO[locusTag]

        for x in elems:
            classCounts[x] += 1

    minValue = -1
    if len(classCounts) < 5:
        minValue = -1

    sorted_x = sorted(classCounts.items(), key=operator.itemgetter(1), reverse=True)


    labels = [goID2label[x[0]] + " " + x[0] for x in classCounts.items() if x[1] > minValue]
    values = [x[1] for x in classCounts.items() if x[1] > minValue]

    labels = [goID2label[x[0]] + " " + x[0] for x in sorted_x]
    values = [x[1] for x in sorted_x]

    indexes = np.arange(len(labels))
    width = 1

    fig, ax = plt.subplots(figsize=(20, 10), dpi=80)
    ax.bar(indexes, values, width)

    ax.set_xticks([i for i in range(0, len(labels))])
    ax.set_xticklabels(labels)

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)

    ax.set_title("GO Plot for diff " + str(diff))
    ax.set_xlabel("GO Classes")
    ax.set_ylabel("Occurrences of GO Classes")

    fig.subplots_adjust(bottom=0.6)
    fig.savefig("/mnt/c/ownCloud/data/hpyloriDB/Wdiff" + str(diff) + ".png", bbox_inches=None)#'tight')

    #fig.tight_layout()
    #plt.show()