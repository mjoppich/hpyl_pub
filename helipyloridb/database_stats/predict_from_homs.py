import random
import sys,os

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from porestat.utils.DataFrame import DataFrame

import numpy as np

from sklearn import tree, svm


def printModelSels(treeModel):
    importances = treeModel.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    showFeatCount = min([30, len(indices), len(homClusterIDs)])

    interestHomCluster = []

    for f in range(showFeatCount):
        interestHomCluster.append(homClusterIDs[f])

        print("%d. feature %d (%f) %s" % (f + 1, indices[f], importances[indices[f]], homClusterIDs[f]))

    for homID in interestHomCluster:
        val = homDB.get_homology_cluster(homID)

        mcc = sum([1 for x in val if x in mc])
        nmcc = sum([1 for x in val if x in nmc])

        print(homID, "MC", mcc, "NMC", nmcc)

    print([int(x.replace('HOMID', '')) for x in interestHomCluster])


if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/combed")
    genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)

    """
    for combid in homDB.combinations:
        elems = homDB.combinations[combid]

        homDB.homologies[combid] = elems

    homDB.finalize()
    homDB.save_to_file(fileLocation + "combed")
    """

    for orgname in homDB.get_all_organisms():
        genomDB.loadGenome(orgname)
    allorgs = list(homDB.get_all_organisms())

    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc] # and not x.startswith("15")


    allData = DataFrame()


    allData.addColumns(allorgs)
    homClusterIDs = []


    for homid in homDB.homologies:

        val = homDB.get_homology_cluster(homid)

        maxlength = 0
        for org in val:

            geneid = val[org]
            seq = genomDB.get_sequence(org, geneid)

            if len(seq) > maxlength:
                maxlength = len(seq)

        if maxlength < 80:
            continue

        allincluded = 0

        if len(val) < 2:
            continue

        for org in nmc:
            if org in val:
                allincluded += 1

        if allincluded > 3:# len(allorgs):
            continue

        allincluded = 0

        for org in mc:
            if org in val:
                allincluded += 1

        if allincluded <= 2:
            continue

        homClusterIDs.append(homid)

    print(allorgs)
    print(len(homClusterIDs))

    orgMatrix = {}

    for org in mc + nmc:

        orgRes = []
        for homID in homClusterIDs:
            val = org in homDB.get_homology_cluster(homID)

            if val:
                orgRes.append(1)
            else:
                orgRes.append(0)

        orgMatrix[org] = orgRes



    trainMatrix = []
    trainRes = []
    trainOrgs = []

    trainSamples = 3

    trainidx = []


    while len(trainOrgs) < 2*trainSamples:

        mcorg = random.choice(mc)
        nmcorg = random.choice(nmc)

        if mcorg in trainOrgs or nmcorg in trainOrgs:
            continue

        mcorgidx = (mc+nmc).index(mcorg)
        nmcorgidx = (mc+nmc).index(nmcorg)

        trainidx += [mcorgidx, nmcorgidx]


        trainMatrix += [orgMatrix[mcorg], orgMatrix[nmcorg]]
        trainRes += [1, 0]
        trainOrgs += [mcorg, nmcorg]

    print(trainOrgs)

    Dforest = np.matrix([orgMatrix[x] for x in mc+nmc])
    colors = [1 for x in mc] + [2 for x in nmc]

    # Build a forest and compute the feature importances
    forest = ExtraTreesClassifier(n_estimators=20, max_features=10, random_state=0)

    forest.fit(Dforest, colors)
    importances = forest.feature_importances_

    model = SelectFromModel(forest, prefit=True)
    X_new = model.transform(Dforest)

    printModelSels(forest)

    print(X_new.shape)

    allrows = []
    for idx in trainidx:
        allrows.append(X_new[idx,])

    D = np.array(allrows)

    clf = tree.DecisionTreeClassifier()
    clf = svm.NuSVC()
    clf = clf.fit(D, trainRes)

    for org in mc+nmc:

        if org in trainOrgs:
            continue

        orgidx = (mc+nmc).index(org)


        predClass = forest.predict(X_new[orgidx,])
        print(org, predClass, org in mc)