import sys
import os
from collections import Counter, defaultdict

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import random
import string
import math

from analysis.homologybuilder import HomologyBuilder
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
import numpy

if __name__ == '__main__':


    homolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    genomeDB = GenomeDB(fileLocation + "/genomes/")

    allowedOrgs = ['CP001217', 'AE000511']

    compareAA = (['W'], ['F', 'G', 'A'])
    compareAA = (['W'], ['H', 'F', 'Y', 'P', 'K'])
    #compareAA = (['W', 'M'], ['H', 'F', 'Y', 'P', 'K'])
    #compareAA = (['W', 'M'], ['F', 'G', 'A'])




    def calculateDifferences(orgI, orgJ, allAA):

        allDiffs = list()
        foundGenes = 0

        printLocusTags = defaultdict(set)

        for homID in homolDB.homologies:

            homCluster = homolDB.homologies[homID]

            orgIT = [x for x in homCluster if x[0] == orgI]
            orgJT = [x for x in homCluster if x[0] == orgJ]

            if len(orgIT) != 1 or len(orgJT) != 1:
                continue

            foundGenes += 1

            orgIT = orgIT[0]
            orgJT = orgJT[0]

            seqI = genomeDB.get_sequence( orgIT[0], orgIT[1] )
            seqJ = genomeDB.get_sequence( orgJT[0], orgJT[1] )


            for diffAA in allAA:
                aaCountI = seqI.count(diffAA)
                aaCountJ = seqJ.count(diffAA)

                countDiff = aaCountI-aaCountJ
                lengthDiff = len(seqI)- len(seqJ)

                kdaDiff = lengthDiff * 0.11

                kdaAvg = (len(seqI)+len(seqJ))/2.0
                kdaAvg *= 0.11

                if abs(kdaDiff) > 5 or (countDiff == 0):
                    continue

                pseudoCount = 0.1
                logFC = float(math.log((aaCountI+pseudoCount)/(aaCountJ+pseudoCount)))

                if math.isnan(logFC):
                    continue

                if orgI in allowedOrgs and orgJ in allowedOrgs:

                    myCountDiff = int(countDiff)
                    if countDiff < 0:
                        myCountDiff = 0-myCountDiff


                    printLocusTags[myCountDiff].add(orgIT[1])

                allDiffs.append( (countDiff, kdaAvg, diffAA) )

        print(orgI, orgJ, foundGenes, len(allDiffs))


        if len(printLocusTags) > 0 and False:

            for absDiff in printLocusTags:

                print(allAA, absDiff, len(printLocusTags[absDiff]), sorted(printLocusTags[absDiff]))

        return allDiffs



    print("Available Orgs:", homolDB.get_all_organisms())
    allOrgs = sorted([x for x in homolDB.get_all_organisms() if x in allowedOrgs])

    for org in allOrgs:
        if org in allowedOrgs:
            genomeDB.loadGenome(org)

    testDiffs = {}
    controlDiffs = {}

    for i in range(0, len(allOrgs)):
        for j in range(i+1, len(allOrgs)):

            orgI = allOrgs[i]
            orgJ = allOrgs[j]

            if not orgI in allowedOrgs or not orgJ in allowedOrgs:
                continue

            allOrgDiffs = calculateDifferences(orgI, orgJ, compareAA[0])
            testDiffs[(orgI, orgJ)] = allOrgDiffs

            allOrgDiffs = calculateDifferences(orgI, orgJ, compareAA[1])
            controlDiffs[(orgI, orgJ)] = allOrgDiffs



    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    import matplotlib.cm as cm
    import numpy as np


    def plotDiffs( diffs, title ):

        colors = iter(cm.rainbow(np.linspace(0, 1, len(diffs))))
        handles = []

        for orgs in diffs:

            orgData = diffs[orgs]

            allX = []
            allY = []

            for diff in orgData:

                allX.append(diff[1])
                allY.append(diff[0])

            h = ax.scatter(allX, allY, color=next(colors), label=str(orgs) + " n=" + str(len(orgData)))
            handles.append(h)


        ax.legend()
        ax.set_title("Scatter plot of kda by count AA for " + title)
        ax.set_xlabel('Protein Difference (kDA)')
        ax.set_ylabel('AA count difference')

        plt.show()

    labelTest = str(compareAA[0])
    labelControl = str(compareAA[1])


    plotDiffs(testDiffs, "TestDiffs: " + labelTest)
    plotDiffs(controlDiffs, "ControlDiffs: " + labelControl)




    from scipy.stats import ks_2samp


    for orgs in testDiffs:

        fig, ax = plt.subplots()
        colors = iter(cm.rainbow(np.linspace(0, 1, 2)))

        orgData = testDiffs[orgs]
        controlData = controlDiffs[orgs]

        ksTest = ks_2samp([x[0] for x in orgData], [x[0] for x in controlData])
        print(orgs, len(orgData), len(controlData))
        print(ksTest)

        ax.hist( [x[0] for x in orgData], 100, histtype='step', cumulative=1, normed = 1, stacked = False, label=str(orgs) + " n=" + str(len(orgData)) + " " + labelTest, color=next(colors) )
        ax.hist( [x[0] for x in controlData], 100, histtype='step', cumulative=1, normed = 1, stacked = False, label=str(orgs) + " n=" + str(len(controlData)) + " " + labelControl, color=next(colors) )


        ax.legend()
        plt.show()