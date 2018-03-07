import sys
import os
from collections import Counter, defaultdict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from scipy.stats import ks_2samp

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random
import string
import math

from analysis.homologybuilder import HomologyBuilder
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
import numpy

interestComps = {(('AE000511', 'HP_0868'), ('CP001217', 'HPP12_0868')),
                 (('AE000511', 'HP_0036'), ('CP001217', 'HPP12_0032')),
                 (('AE000511', 'HP_0963'), ('CP001217', 'HPP12_0958')),
                 (('AE000511', 'HP_1282'), ('CP001217', 'HPP12_1248')),
                 (('AE000511', 'HP_0568'), ('CP001217', 'HPP12_0574')),
                 (('AE000511', 'HP_0286'), ('CP001217', 'HPP12_0285')),
                 (('AE000511', 'HP_0519'), ('CP001217', 'HPP12_0525')),
                 (('AE000511', 'HP_0104'), ('CP001217', 'HPP12_0106')),
                 (('AE000511', 'HP_0091'), ('CP001217', 'HPP12_0094')),
                 (('AE000511', 'HP_0342'), ('CP001217', 'HPP12_0337')),
                 (('AE000511', 'HP_0656'), ('CP001217', 'HPP12_0669')),
                 (('AE000511', 'HP_0043'), ('CP001217', 'HPP12_0038')),
                 (('AE000511', 'HP_0108'), ('CP001217', 'HPP12_0110')),
                 (('AE000511', 'HP_1213'), ('CP001217', 'HPP12_1179')),
                 (('AE000511', 'HP_1105'), ('CP001217', 'HPP12_1070')),
                 (('AE000511', 'HP_0661'), ('CP001217', 'HPP12_0674')),
                 (('AE000511', 'HP_0430'), ('CP001217', 'HPP12_0992')),
                 (('AE000511', 'HP_0048'), ('CP001217', 'HPP12_0042')),
                 (('AE000511', 'HP_0860'), ('CP001217', 'HPP12_0860')),
                 (('AE000511', 'HP_0558'), ('CP001217', 'HPP12_0565')),
                 (('AE000511', 'HP_0289'), ('CP001217', 'HPP12_0288')),
                 (('AE000511', 'HP_0743'), ('CP001217', 'HPP12_0752')),
                 (('AE000511', 'HP_0744_2'), ('CP001217', 'HPP12_0753')),
                 (('AE000511', 'HP_0312'), ('CP001217', 'HPP12_0311')),
                 (('AE000511', 'HP_1046'), ('CP001217', 'HPP12_0398')),
                 (('AE000511', 'HP_1184'), ('CP001217', 'HPP12_1149')),
                 (('AE000511', 'HP_0755'), ('CP001217', 'HPP12_0765')),
                 (('AE000511', 'HP_0818'), ('CP001217', 'HPP12_0825')),
                 (('AE000511', 'HP_1533'), ('CP001217', 'HPP12_1507')),
                 (('AE000511', 'HP_0235'), ('CP001217', 'HPP12_0235')),
                 (('AE000511', 'HP_1141'), ('CP001217', 'HPP12_1107')),
                 (('AE000511', 'HP_0099'), ('CP001217', 'HPP12_0101')),
                 (('AE000511', 'HP_1455'), ('CP001217', 'HPP12_1433')),
                 (('AE000511', 'HP_0566'), ('CP001217', 'HPP12_0572')),
                 (('AE000511', 'HP_0759'), ('CP001217', 'HPP12_0769')),
                 (('AE000511', 'HP_0783'), ('CP001217', 'HPP12_0792')),
                 (('AE000511', 'HP_0897'), ('CP001217', 'HPP12_0892')),
                 (('AE000511', 'HP_1277'), ('CP001217', 'HPP12_1243')),
                 (('AE000511', 'HP_0599'), ('CP001217', 'HPP12_0606')),
                 (('AE000511', 'HP_0064'), ('CP001217', 'HPP12_0065')),
                 (('AE000511', 'HP_1349'), ('CP001217', 'HPP12_1312')),
                 (('AE000511', 'HP_0142'), ('CP001217', 'HPP12_0141')),
                 (('AE000511', 'HP_1229'), ('CP001217', 'HPP12_1195')),
                 (('AE000511', 'HP_0270'), ('CP001217', 'HPP12_0269')),
                 (('AE000511', 'HP_1107'), ('CP001217', 'HPP12_1072')),
                 (('AE000511', 'HP_0793'), ('CP001217', 'HPP12_0800')),
                 (('AE000511', 'HP_1072'), ('CP001217', 'HPP12_0372')),
                 (('AE000511', 'HP_0415'), ('CP001217', 'HPP12_1009')),
                 (('AE000511', 'HP_1100'), ('CP001217', 'HPP12_1065')),
                 (('AE000511', 'HP_0965'), ('CP001217', 'HPP12_0960')),
                 (('AE000511', 'HP_0762'), ('CP001217', 'HPP12_0772'))}

print(len(interestComps))

if __name__ == '__main__':

    observeChanges = [('W', 'L'), ('W', 'C')]


    def calculateDifferences(homolDB, orgI, orgJ,):

        obsChanges = {'interest': defaultdict(list), 'others': defaultdict(list)}
        listenToChanges = [x[0] for x in observeChanges]

        orgSubMatrix = Counter()
        orgSubMatRel = Counter()

        aaLength = 0

        for homID in homolDB.homologies:

            homCluster = homolDB.homologies[homID]


            orgIT = [x for x in homCluster if x[0] == orgI]
            orgJT = [x for x in homCluster if x[0] == orgJ]

            if len(orgIT) != 1 or len(orgJT) != 1:
                continue

            orgIT = orgIT[0]
            orgJT = orgJT[0]


            seqI = genomeDB.get_sequence( orgIT[0], orgIT[1] )
            seqJ = genomeDB.get_sequence( orgJT[0], orgJT[1] )

            if abs(len(seqI)-len(seqJ)) > 15:
                continue

            if '*' in seqI+seqJ:
                continue

            if 'X' in seqI+seqJ:
                continue

            if 'U' in seqI+seqJ:
                continue

            #print("HOMID", orgI, orgJ, homID)


            aaLength += (len(seqI)+len(seqJ))/2.0
            alignments = pairwise2.align.globalds(seqI, seqJ, matrix, -50, -.2)

            alignment = alignments[0]

            alignI = alignment[0]
            alignJ = alignment[1]

            printAlignment=False

            obsCount = Counter()
            refCount = 0
            refCountChanges = Counter()

            for i in range(0, len(alignI)):

                aaI = alignI[i]
                aaJ = alignJ[i]

                if aaI in listenToChanges or aaJ in listenToChanges:
                    refCount += 1

                if aaI in listenToChanges and aaI == aaJ:
                    refCountChanges[aaI] += 1

                if (aaI, aaJ) in observeChanges:
                    obsCount[(aaI, aaJ)] += 1

                if (aaJ, aaI) in observeChanges:
                    obsCount[(aaJ, aaI)] += 1


            if printAlignment:
                print("WL Changes", orgIT, orgJT, orgSubMatrix[('L', 'W')])
                print(alignment[0])
                print(alignment[1])

            accept = refCount > 0

            if len(obsCount) > 0:
                print(orgIT, orgJT, obsCount, refCount, refCountChanges, accept, (orgIT, orgJT) in interestComps or (orgJT, orgIT) in interestComps)

            if not accept:
                continue

            for chg in observeChanges:
                frac = obsCount[chg] / aaLength
                fracRef = obsCount[chg] / refCount if refCount > 0 else 0


                if (orgIT, orgJT) in interestComps or (orgJT, orgIT) in interestComps:
                    obsChanges['interest'][chg].append( (frac, fracRef) )
                else:
                    obsChanges['others'][chg].append( (frac, fracRef) )

        return obsChanges


    hpHomolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    cbHomolDB = HomologyDatabase.loadFromFile(fileLocation + "../cbdb/"+ "/cbj")


    genomeDB = GenomeDB(fileLocation + "/genomes/")


    genomeDB.loadGenome(fileLocation + "/genomes/CP001217.gb")
    genomeDB.loadGenome(fileLocation + "/genomes/AE000511.gb")

    genomeDB.fileExtension = '.gbff'
    genomeDB.fileFormat = 'gb'

    genomeDB.loadGenome(fileLocation + "../cbdb/genomes/NC003912.gbff")
    genomeDB.loadGenome(fileLocation + "../cbdb/genomes/NC002163.gbff")

    matrix = matlist.blosum80

    subMatrix = {}


    print("Starting HPP")
    (orgI, orgJ) = ('AE000511', 'CP001217')

    observedChanges = calculateDifferences(hpHomolDB, orgI, orgJ)

    for chg in observeChanges:

        cnt = 0

        for group in observedChanges:

            for elem in observedChanges[group][chg]:
                if elem[0] > 0.0:
                    cnt += 1

        print(chg, cnt)


    def analyseChanges( observedChanges):

        for obschg in observeChanges:
            interestData = [x[1] for x in observedChanges['interest'][obschg]]
            genomeData = [x[1] for x in observedChanges['others'][obschg]]

            print(interestData)
            print(genomeData)

            ksTestInt = ks_2samp(interestData, genomeData)

            print(obschg, ksTestInt)

            fig, ax = plt.subplots()
            colors = iter(cm.rainbow(np.linspace(0, 1, 3)))

            ax.hist(interestData, 100, histtype='step', cumulative=1, normed=1, stacked=False,
                    label=" n=" + str(len(interestData)) + " Interst", color=next(colors))
            ax.hist(genomeData, 100, histtype='step', cumulative=1, normed=1, stacked=False,
                    label=" n=" + str(len(genomeData)) + " Others", color=next(colors))

            ax.set_title("Distribution of " + str(obschg) + " " + str(ksTestInt.pvalue))

            ax.legend()
            plt.show()


    analyseChanges(observedChanges)

    exit(0)

    print("Starting Differences")
    (orgI, orgJ) = ('NC003912', 'NC002163')

    observedChanges = calculateDifferences(cbHomolDB, orgI, orgJ)

    analyseChanges(observedChanges)
