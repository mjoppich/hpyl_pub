import sys
import os
from collections import Counter, defaultdict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils.utils import fileLocation

if __name__ == '__main__':


    homolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    genomeDB = GenomeDB(fileLocation + "/genomes/")

    allowedOrgs = ['CP001217', 'AE000511']

    compareAA = (['W'], ['F', 'G', 'A'])
    compareAA = (['W'], ['F', 'Y'])
    #compareAA = (['W', 'M'], ['H', 'F', 'Y', 'P', 'K'])
    #compareAA = (['W', 'M'], ['F', 'G', 'A'])
    #compareAA = (['W'], ['F', 'Y', 'W' ])
    #compareAA = (['W'], ['F', 'Y'])



    aaUsage = {}

    aaUsage['CP001217'] = {'A': 0.06820758269610876, 'W': 0.007104265486069513, 'E': 0.06884757145728172,
                           'N': 0.059195058038005184, 'I': 0.07177239814325212, 'C': 0.010989075313799368,
                           'M': 0.02266106546421624, 'T': 0.043905570438763024, 'S': 0.06825245995680077,
                           'Y': 0.03682667036091074, 'K': 0.08927257862788751, 'H': 0.020990850892374572,
                           'V': 0.05570438763026601, 'Q': 0.036988618736451456, 'L': 0.11184779194121469,
                           'G': 0.05795995777635125, 'F': 0.05385856638615127, 'D': 0.047431362180956824,
                           'P': 0.03316624683664092, 'R': 0.03470768318214899, '*': 0.0003102384543490944}

    aaUsage['AE000511'] = {'E': 0.06852955118080811, 'I': 0.07175324483746802, 'M': 0.02275664161534785,
                           'T': 0.04388275117372739, 'S': 0.0684587439985996, 'K': 0.08920524838569459,
                           'H': 0.0212185522684851, 'L': 0.11177897144700377, 'P': 0.033005981240030444,
                           'F': 0.05397080777226837, 'D': 0.04761979690139903, 'R': 0.03493744382138495,
                           'W': 0.0071593928677498855, 'A': 0.06828369290925076, 'N': 0.05843559398374975,
                           '*': 0.0002006203495907935, 'C': 0.011122628205254287, 'X': 3.7370457276716434e-05,
                           'V': 0.0558530986993114, 'Q': 0.037242610975506615, 'G': 0.05770982036611247,
                           'Y': 0.03683743654398011}

    aaUsage['CP001217'] = Counter({'L': 0.11358863299560508, 'K': 0.0874111952781647, 'I': 0.07241631526391779, 'A': 0.07042556303942708, 'E': 0.06859127249136603, 'S': 0.06804047650848773, 'G': 0.05978235290376336, 'V': 0.0583907344158446, 'N': 0.05513302422611892, 'F': 0.053709604584454326, 'D': 0.04633428100771496, 'T': 0.04311600424863415, 'Y': 0.03615536771674076, 'Q': 0.03585898096382937, 'R': 0.03454750138335021, 'P': 0.03356293766337842, 'M': 0.023509956941237826, 'H': 0.021657857747079068, 'C': 0.010856913888835888, 'W': 0.00691102673204984, '-': 0.0015989620103417355})
    aaUsage['AE000511'] = aaUsage['CP001217']


    #compareAA = (compareAA[0], [x for x in aaUsage['CP001217'] if x != '*' and not x in compareAA[0]])


    trpChange = ([], [])

    targetLT = []

    def calculateDifferences(orgI, orgJ, allAA, only_diffs=False, excludeGenes=set()):

        allDiffs = list()
        foundGenes = 0
        genesWithDiff = set()

        genePairWithAA = set()

        printLocusTags = defaultdict(set)

        trpCounts = Counter()

        for homID in homolDB.homologies:

            homCluster = homolDB.homologies[homID]

            orgIT = [x for x in homCluster if x[0] == orgI]
            orgJT = [x for x in homCluster if x[0] == orgJ]

            if len(orgIT) != 1 or len(orgJT) != 1:
                continue

            foundGenes += 1

            orgIT = orgIT[0]
            orgJT = orgJT[0]

            if orgIT[1] in excludeGenes:
                continue

            seqI = genomeDB.get_sequence( orgIT[0], orgIT[1] )
            seqJ = genomeDB.get_sequence( orgJT[0], orgJT[1] )

            if '*' in seqI or '*' in seqJ:
                continue

            if 'X' in seqI or 'X' in seqJ:
                continue

            for aa in allAA:
                if aa in seqI or aa in seqJ:
                    genePairWithAA.add(orgIT[1])

            trpCounts[orgIT[0]] += seqI.count('W')
            trpCounts[orgJT[0]] += seqJ.count('W')


            for diffAA in allAA:

                diffGenes = set([x[3] for x in allDiffs])

                if orgIT[1] in diffGenes:
                    continue

                aaCountI = seqI.count(diffAA)
                aaCountJ = seqJ.count(diffAA)

                countDiff = aaCountI-aaCountJ
                lengthDiff = len(seqI)- len(seqJ)

                kdaDiff = lengthDiff * 0.11

                kdaAvg = (len(seqI)+len(seqJ))/2.0
                kdaAvg *= 0.11

                if abs(kdaDiff) > 1000 or (countDiff == 0):
                    continue

                if not (aaCountI == 0 or aaCountJ == 0):
                    continue

                genesWithDiff.add(orgIT[1])

                pseudoCount = 0.1
                logFC = float(math.log((aaCountI+pseudoCount)/(aaCountJ+pseudoCount)))

                if math.isnan(logFC):
                    continue

                add_element = False

                if not only_diffs and not add_element:

                    matrix = matlist.blosum80
                    alignments = pairwise2.align.globalds(seqI, seqJ, matrix, -50, -.2)

                    alignment = alignments[0]

                    changes = 0
                    losses = 0
                    changesAAA = 0

                    printAlignment = False

                    for i in range(0, len(alignment[0])):

                        if not alignment[0][i] == alignment[1][i]:

                            aaI = alignment[0][i]
                            aaJ = alignment[0][j]

                            if set([aaI, aaJ]) == set(['W', 'L']):
                                printAlignment = True

                            if alignment[0][i] == diffAA and alignment[1][i] != '-':
                                changesAAA += 1

                            if alignment[1][i] == '-':
                                losses += 1
                            else:
                                changes += 1



                    if printAlignment and False:
                        print("WL Changes", orgIT)
                        print(alignment[0])
                        print(alignment[1])


                    expectedWChange = changes * aaUsage[orgI][diffAA]
                    expectedWLosses = losses * aaUsage[orgI][diffAA]

                    actualChange = expectedWChange-expectedWLosses

                    expectedCount = (aaCountI-changesAAA) + expectedWChange-expectedWLosses

                    expCountDiff = abs(expectedCount-aaCountJ)

                    if (abs(actualChange) >= 0.2 and abs(expCountDiff) > 1) or (abs(actualChange) < 0.2 and abs(expCountDiff) > 0.8):
                        #print(diffAA, orgIT, len(seqI), aaCountI, orgJT, len(seqJ), aaCountJ, kdaDiff, changes, changesAAA, countDiff, expectedWChange, expectedCount, expCountDiff)
                        targetLT.append(orgIT[1])
                        add_element = True
                    else:
                        pass#print(diffAA, orgIT, len(seqI), aaCountI, orgJT, len(seqJ), aaCountJ, kdaDiff, changes, changesAAA, countDiff, expectedWChange, expectedCount, expCountDiff)


                    trpChange[0].append(countDiff)
                    trpChange[1].append(changes*aaUsage[orgI][diffAA])

                else:
                    add_element = True


                if add_element:
                    myCountDiff = int(countDiff)
                    if countDiff < 0:
                        myCountDiff = 0-myCountDiff

                    printLocusTags[myCountDiff].add(orgIT[1])

                    allDiffs.append( (countDiff, kdaAvg, diffAA, orgIT, orgJT) )

        print(orgI, orgJ, foundGenes, len(allDiffs))


        if len(printLocusTags) > 0 and False:
            for absDiff in printLocusTags:
                print(allAA, absDiff, len(printLocusTags[absDiff]), sorted(printLocusTags[absDiff]))

        diffGenes = set([x[3][1] for x in allDiffs])
        diffGenePairs = set([(x[3], x[4]) for x in allDiffs])

        print(trpCounts)

        return allDiffs, diffGenes, genesWithDiff, genePairWithAA, diffGenePairs


    print("Available Orgs:", homolDB.get_all_organisms())
    allOrgs = sorted([x for x in homolDB.get_all_organisms() if x in allowedOrgs])
    print("Used Orgs:", allOrgs)


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

            allOrgDiffs, allDiffGenes, genesWithDiff, genePairsWithAA, diffGenePairs = calculateDifferences(orgI, orgJ, compareAA[0])
            testDiffs[(orgI, orgJ)] = allOrgDiffs

            print("DiffGenes:", len(allDiffGenes))
            print("GenesWithDiff:", len(genesWithDiff))
            print("GenePairsWithAA:", len(genePairsWithAA))
            print("DiffGenes:", allDiffGenes)
            print("GenesWithDiff", genesWithDiff)
            print("DiffGenePairs", diffGenePairs)

            allOrgDiffs, allDiffGenes, genesWithDiff, genePairsWithAA, diffGenePairs = calculateDifferences(orgI, orgJ, compareAA[1], only_diffs=False, excludeGenes=allDiffGenes)
            controlDiffs[(orgI, orgJ)] = allOrgDiffs
            print("DiffGenes:", len(allDiffGenes))
            print("GenesWithDiff:", len(genesWithDiff))
            print("GenePairsWithAA:", len(genePairsWithAA))
            print("DiffGenePairs", diffGenePairs)


            fig, ax = plt.subplots()
            ax.scatter(trpChange[0], trpChange[1])
            plt.show()

    exit(0)


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