from collections import defaultdict, Counter
from itertools import chain

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils.utils import fileLocation

import Bio.Alphabet.Reduced


def countAAFractions(homDB, orgI, orgJ):

    aaAbsCounts = defaultdict(lambda: defaultdict(list))
    aaRelCounts = defaultdict(lambda: defaultdict(list))

    ignoredSeqPairs = 0

    for homID in homDB.homologies:
        homGroup = homDB.get_homology_cluster(homID)
        homCluster = [(x, homGroup[x]) for x in homGroup]

        orgIT = [x for x in homCluster if x[0] == orgI]
        orgJT = [x for x in homCluster if x[0] == orgJ]

        if len(orgIT) == 0 or len(orgJT) == 0:
            continue

        orgIT = orgIT[0]
        orgJT = orgJT[0]

        seqI = genomeDB.get_sequence(orgIT[0], orgIT[1])
        seqJ = genomeDB.get_sequence(orgJT[0], orgJT[1])

        if '*' in seqI or '*' in seqJ:
            ignoredSeqPairs += 1
            continue

        if 'X' in seqI or 'X' in seqJ:
            ignoredSeqPairs += 1
            continue

        if 'U' in seqI or 'U' in seqJ:
            ignoredSeqPairs += 1
            continue

        aaSeqCounterI = Counter(seqI)
        aaSeqCounterJ = Counter(seqJ)

        for aa in allAA:
            aaAbsCounts[orgIT[0]][aa].append(aaSeqCounterI[aa])
            aaRelCounts[orgIT[0]][aa].append(aaSeqCounterI[aa]/len(seqI))

        for aa in allAA:
            aaAbsCounts[orgJT[0]][aa].append(aaSeqCounterJ[aa])
            aaRelCounts[orgJT[0]][aa].append(aaSeqCounterJ[aa]/len(seqJ))


    print("Ignored", ignoredSeqPairs)
    return aaAbsCounts, aaRelCounts

if __name__ == '__main__':



    allAA = [x for x in Bio.Alphabet.Reduced.hp_model_tab]
    print(allAA)


    hpHomolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    cbHomolDB = HomologyDatabase.loadFromFile(fileLocation + "../cbdb/"+ "/cbj")


    genomeDB = GenomeDB(fileLocation + "/genomes/")


    genomeDB.loadGenome(fileLocation + "/genomes/CP001217.gb")
    genomeDB.loadGenome(fileLocation + "/genomes/AE000511.gb")

    genomeDB.fileExtension = '.gbff'
    genomeDB.fileFormat = 'gb'

    genomeDB.loadGenome(fileLocation + "../cbdb/genomes/NC003912.gbff")
    genomeDB.loadGenome(fileLocation + "../cbdb/genomes/NC002163.gbff")


    print("Starting HPP")
    (orgI, orgJ) = ('CP001217', 'AE000511')

    absCount, relCount = countAAFractions(hpHomolDB, orgI, orgJ)

    print((orgI, orgJ))

    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import ks_2samp
    import matplotlib.pyplot as plt

    import matplotlib.cm as cm
    import numpy as np


    def plotDistributions( analyseAA = ['W']):
        fig, ax = plt.subplots()
        colors = iter(cm.rainbow(np.linspace(0, 1, 3)))

        genome = str((orgI, orgJ))

        org1data = [relCount[orgI][aa] for aa in analyseAA]
        org2data = [relCount[orgJ][aa] for aa in analyseAA]

        org1data = list(chain.from_iterable(org1data))
        org2data = list(chain.from_iterable(org2data))

        print(len(org1data), len(org2data))



        ksTestInt = ks_2samp(org1data, org2data)

        ax.hist( org1data, 100, histtype='step', cumulative=1, normed = 1, stacked = False, label=str(genome) + " n=" + str(len(org1data)) + " Interst", color=next(colors) )
        ax.hist( org2data, 100, histtype='step', cumulative=1, normed = 1, stacked = False, label=str(genome) + " n=" + str(len(org2data)) + " Others", color=next(colors) )

        ax.set_title(genome + " Distribution of " + str(analyseAA) + " KS-pVal=" + str(ksTestInt.pvalue))

        ax.legend()
        plt.show()



    plotDistributions(['W'])

    print("Starting Differences")
    (orgI, orgJ) = ('NC003912', 'NC002163')

    absCount, relCount = countAAFractions(cbHomolDB, orgI, orgJ)

    print((orgI, orgJ))

    plotDistributions(['W'])
