import random
import re
import sys,os

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.svm import LinearSVC

from database.homDBAnalyser import HomDBAnalyser

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
from porestat.utils.DataFrame import DataFrame, DataRow

import pandas as pd
import seaborn as sns
import numpy as np

from sklearn import decomposition, tree, svm
from sklearn import datasets

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def printModelSels(treeModel):


    importances = treeModel.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    print(importances)

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

    intHomIDs = [x for x in interestHomCluster]
    print(intHomIDs)

    return intHomIDs


if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp_split")
    genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)

    """
    for combid in homDB.combinations:
        elems = homDB.combinations[combid]

        homDB.homologies[combid] = elems

    homDB.finalize()
    homDB.save_to_file(fileLocation + "combed")
    """

    thoms = ['HOMID1448', 'HOMID1742', 'HOMID1692', 'HOMID1795', 'HOMID2024', 'HOMID2027', 'HOMID1621', 'HOMID1338', 'HOMID1672', 'HOMID1693']
    thoms = [1549, 1753, 1539, 1971, 1585, 1547, 1820, 1545, 1544, 1555, 1787, 1776, 1634]


    thoms = ['HOMID2364', 'HOMID2388', 'HOMID1942', 'HOMID2293', 'HOMID2311', 'HOMID2294', 'HOMID2363', 'HOMID2258', 'HOMID2260', 'HOMID2287', 'HOMID1924', 'HOMID1937', 'HOMID1927', 'HOMID2263', 'HOMID2323', 'HOMID2290', 'HOMID2354', 'HOMID1934', 'HOMID2306', 'HOMID2349', 'HOMID2316', 'HOMID1919', 'HOMID2320', 'HOMID2308', 'HOMID1925', 'HOMID2278', 'HOMID2261', 'HOMID2353', 'HOMID2285', 'HOMID2291', 'HOMID2313', 'HOMID2385', 'HOMID2376', 'HOMID2297', 'HOMID2347', 'HOMID2345', 'HOMID2370', 'HOMID2340', 'HOMID2283', 'HOMID2299', 'HOMID2372', 'HOMID1928', 'HOMID2281', 'HOMID2375', 'HOMID2366', 'HOMID2277', 'HOMID1916', 'HOMID2284', 'HOMID1922', 'HOMID2307', 'HOMID1943', 'HOMID2262', 'HOMID2300', 'HOMID2322', 'HOMID2050', 'HOMID2346', 'HOMID2302', 'HOMID2296', 'HOMID2251', 'HOMID2391', 'HOMID2298', 'HOMID2275', 'HOMID2374', 'HOMID2337', 'HOMID2315', 'HOMID2371', 'HOMID2348', 'HOMID2256', 'HOMID2355', 'HOMID2387', 'HOMID2317', 'HOMID2368', 'HOMID2289', 'HOMID2377', 'HOMID1923', 'HOMID1939', 'HOMID2373', 'HOMID2301', 'HOMID2324', 'HOMID2386', 'HOMID2310', 'HOMID2360', 'HOMID2352', 'HOMID2344', 'HOMID2312', 'HOMID2042', 'HOMID2381', 'HOMID2357', 'HOMID1920', 'HOMID2380', 'HOMID2367', 'HOMID2369', 'HOMID1935', 'HOMID2341', 'HOMID2274', 'HOMID2305', 'HOMID2384', 'HOMID2319', 'HOMID2359', 'HOMID2379', 'HOMID2303', 'HOMID2336', 'HOMID2342', 'HOMID2390', 'HOMID2365', 'HOMID2358', 'HOMID2265', 'HOMID2250', 'HOMID2350', 'HOMID2314', 'HOMID2351', 'HOMID2318', 'HOMID2304', 'HOMID1917', 'HOMID2288', 'HOMID2282', 'HOMID2383', 'HOMID1926', 'HOMID2339', 'HOMID2389', 'HOMID2378', 'HOMID1930', 'HOMID2321', 'HOMID2362', 'HOMID1940', 'HOMID2279', 'HOMID2309', 'HOMID2338', 'HOMID2286', 'HOMID2382', 'HOMID1933', 'HOMID2343', 'HOMID2361', 'HOMID2272', 'HOMID2292', 'HOMID2356', 'sp_HOMID2424', 'sp_HOMID2436']

    thoms = ['HOMID403', 'HOMID1649', 'HOMID981', 'HOMID607', 'HOMID346']

    thoms =['HOMID1610', 'HOMID403', 'HOMID1508', 'HOMID1243', 'HOMID1649', 'HOMID1849', 'HOMID1254', 'HOMID1229', 'HOMID1826', 'HOMID1795', 'HOMID1703', 'HOMID1662', 'HOMID1623', 'HOMID981', 'HOMID1573', 'HOMID1936', 'HOMID1609', 'HOMID607', 'HOMID1396', 'HOMID346', 'HOMID1611', 'HOMID1958', 'HOMID1618', 'HOMID1781', 'HOMID1898', 'HOMID1749', 'HOMID1463', 'HOMID1608', 'HOMID1335', 'HOMID1830', 'HOMID1671', 'HOMID1531', 'HOMID2022', 'HOMID1083', 'HOMID1920', 'HOMID1935', 'HOMID1600', 'HOMID632', 'HOMID1845', 'HOMID1235', 'HOMID1866', 'HOMID1633', 'HOMID1739', 'HOMID1637', 'HOMID1323', 'HOMID1357', 'HOMID1612', 'HOMID1582', 'HOMID1677', 'HOMID1688', 'HOMID1815', 'sp_HOMID2404', 'sp_HOMID2436', 'sp_HOMID2463', 'sp_HOMID2486', 'sp_HOMID2494', 'sp_HOMID2501', 'sp_HOMID2502']



    targetHOMIDS = []

    for x in thoms:
        if isinstance(x, str):
            targetHOMIDS.append(x)




    expectedFeatureCount = 5

    for orgname in homDB.get_all_organisms():
        genomDB.loadGenome(orgname)
    allorgs = list(homDB.get_all_organisms())

    extra = ['AE001439', 'CP009259']

    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc and not x in extra and not x.startswith("6_")] # and not x.startswith("15")


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


        if not homid in targetHOMIDS:
            continue



        if len(val)  <2:
            continue

        if maxlength < 80:
            continue

        nmcIncluded = 0
        for org in nmc:
            if org in val:
                nmcIncluded += 1

        mcInlcuded = 0
        for org in mc:
            if org in val:
                mcInlcuded += 1

        if mcInlcuded >= 8 and nmcIncluded == 6:
            continue


        #if nmcIncluded > 2 or mcInlcuded < 2:
            #continue

        """
        #if nmcIncluded < mcInlcuded:
        #    continue
        """

        homClusterIDs.append(homid)

    df = DataFrame()
    df.addColumns(['sample'] + [homid for homid in targetHOMIDS])

    for org in mc+nmc:


        rowdict = {
            'sample': org
        }

        for homid in targetHOMIDS:

            homID = homid
            val = homDB.get_cluster(homID)

            if org in val:
                rowdict[homID] = 1
            else:
                rowdict[homID] = 0

        dfrow = DataRow.fromDict(rowdict)
        df.addRow(dfrow)

    df.export(outFile=None)

    print(allorgs)
    print(len(homClusterIDs))

    if len(homClusterIDs) < 10:
        print(homClusterIDs)

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


    numCorrect = 0
    targetHOMS = []

    tries = 0
    while numCorrect < len(mc+nmc) and tries < 50:

        tries += 1

        trainMatrix = []
        trainRes = []
        trainOrgs = []

        trainSamples = 4

        trainidx = []
        numCorrect = 0


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

        Dall = np.matrix([orgMatrix[x] for x in mc+nmc])

        trainRows = []
        for idx in trainidx:
            trainRows.append(orgMatrix[(mc+nmc)[idx]])

        trainingMat = np.array(trainMatrix)

        # Build a forest and compute the feature importances
        forest = ExtraTreesClassifier(n_estimators=5, max_features=expectedFeatureCount)

        #forest = LinearSVC()

        forest.fit(trainingMat, trainRes)

        model = SelectFromModel(forest, prefit=True)
        X_new = model.transform(Dall)

        print(X_new.shape)

        print("Training")
        for org in mc+nmc:

            if not org in trainOrgs:
                continue

            orgidx = (mc+nmc).index(org)

            predClass = forest.predict(Dall[orgidx,])
            predProbs = 0#forest.predict_proba(Dall[orgidx,])

            isCorrect = (org in mc) == (1 in predClass)

            if isCorrect:
                numCorrect += 1

            print(org, predClass, predProbs, org in mc, isCorrect)

        print("Prediction")

        for org in mc+nmc:

            if org in trainOrgs:
                continue

            orgidx = (mc+nmc).index(org)

            predClass = forest.predict(Dall[orgidx,])
            predProbs = 0#forest.predict_proba(Dall[orgidx,])

            isCorrect = (org in mc) == (1 in predClass)

            if isCorrect:
                numCorrect += 1

            print(org, predClass, predProbs, org in mc, isCorrect)




    if isinstance(forest, LinearSVC):

        allhoms = []

        for  (idx, x) in enumerate(forest.coef_[0]):

            if abs(x) > 0.1:
                print(idx, homClusterIDs[idx], x)
                allhoms.append(homClusterIDs[idx])


        print(allhoms)

        targetHOMS = [x for x in allhoms]

        print(targetHOMS)


    else:



        targetHOMS = printModelSels(forest)
        importances = forest.feature_importances_

        analyse = HomDBAnalyser(homDB, genomDB)

        def printHOM(homid):
            print(homid)

            aligned = analyse.cluster_align(homid)
            longest = ""

            for rec in aligned._records:

                if len(str(rec.seq).replace('-', '')) > len(longest):
                    longest = str(rec.seq).replace('-', '')

                print(rec.seq, rec.id)

            return (homid, longest)


        allseqs = []

        for homid in targetHOMS:
            allseqs.append(printHOM(homid))



    df = DataFrame()
    df.addColumns(['sample'] + [homid for homid in targetHOMS])

    for org in mc+nmc:


        rowdict = {
            'sample': org
        }

        for homid in targetHOMS:

            homID = homid
            val = homDB.get_cluster(homID)

            if org in val:
                rowdict[homID] = 1
            else:
                rowdict[homID] = 0

        dfrow = DataRow.fromDict(rowdict)
        df.addRow(dfrow)

    df.export(outFile=None)


    print(forest.get_params())

    for elem in allseqs:
        print(">" + elem[0] + " " + str(len(elem[1])))
        print(elem[1])