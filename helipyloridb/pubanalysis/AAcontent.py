from collections import defaultdict, OrderedDict

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib
import numpy as np


if __name__ == '__main__':

    homolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    genomeDB = GenomeDB(fileLocation + "/genomes/")

    allowedOrgs = ['AE000511', 'CP001217', 'CP001173', 'AE001439']
    allAA = ['W']


    diffGenePairs = {(('AE000511', 'HP_0868'), ('CP001217', 'HPP12_0868')), (('AE000511', 'HP_0036'), ('CP001217', 'HPP12_0032')),
     (('AE000511', 'HP_0963'), ('CP001217', 'HPP12_0958')), (('AE000511', 'HP_1282'), ('CP001217', 'HPP12_1248')),
     (('AE000511', 'HP_0568'), ('CP001217', 'HPP12_0574')), (('AE000511', 'HP_0286'), ('CP001217', 'HPP12_0285')),
     (('AE000511', 'HP_0519'), ('CP001217', 'HPP12_0525')), (('AE000511', 'HP_0104'), ('CP001217', 'HPP12_0106')),
     (('AE000511', 'HP_0091'), ('CP001217', 'HPP12_0094')), (('AE000511', 'HP_0342'), ('CP001217', 'HPP12_0337')),
     (('AE000511', 'HP_0656'), ('CP001217', 'HPP12_0669')), (('AE000511', 'HP_0043'), ('CP001217', 'HPP12_0038')),
     (('AE000511', 'HP_0108'), ('CP001217', 'HPP12_0110')), (('AE000511', 'HP_1213'), ('CP001217', 'HPP12_1179')),
     (('AE000511', 'HP_1105'), ('CP001217', 'HPP12_1070')), (('AE000511', 'HP_0661'), ('CP001217', 'HPP12_0674')),
     (('AE000511', 'HP_0430'), ('CP001217', 'HPP12_0992')), (('AE000511', 'HP_0048'), ('CP001217', 'HPP12_0042')),
     (('AE000511', 'HP_0860'), ('CP001217', 'HPP12_0860')), (('AE000511', 'HP_0558'), ('CP001217', 'HPP12_0565')),
     (('AE000511', 'HP_0289'), ('CP001217', 'HPP12_0288')), (('AE000511', 'HP_0743'), ('CP001217', 'HPP12_0752')),
     (('AE000511', 'HP_0744_2'), ('CP001217', 'HPP12_0753')), (('AE000511', 'HP_0312'), ('CP001217', 'HPP12_0311')),
     (('AE000511', 'HP_1046'), ('CP001217', 'HPP12_0398')), (('AE000511', 'HP_1184'), ('CP001217', 'HPP12_1149')),
     (('AE000511', 'HP_0755'), ('CP001217', 'HPP12_0765')), (('AE000511', 'HP_0818'), ('CP001217', 'HPP12_0825')),
     (('AE000511', 'HP_1533'), ('CP001217', 'HPP12_1507')), (('AE000511', 'HP_0235'), ('CP001217', 'HPP12_0235')),
     (('AE000511', 'HP_1141'), ('CP001217', 'HPP12_1107')), (('AE000511', 'HP_0099'), ('CP001217', 'HPP12_0101')),
     (('AE000511', 'HP_1455'), ('CP001217', 'HPP12_1433')), (('AE000511', 'HP_0566'), ('CP001217', 'HPP12_0572')),
     (('AE000511', 'HP_0759'), ('CP001217', 'HPP12_0769')), (('AE000511', 'HP_0783'), ('CP001217', 'HPP12_0792')),
     (('AE000511', 'HP_0897'), ('CP001217', 'HPP12_0892')), (('AE000511', 'HP_1277'), ('CP001217', 'HPP12_1243')),
     (('AE000511', 'HP_0599'), ('CP001217', 'HPP12_0606')), (('AE000511', 'HP_0064'), ('CP001217', 'HPP12_0065')),
     (('AE000511', 'HP_1349'), ('CP001217', 'HPP12_1312')), (('AE000511', 'HP_0142'), ('CP001217', 'HPP12_0141')),
     (('AE000511', 'HP_1229'), ('CP001217', 'HPP12_1195')), (('AE000511', 'HP_0270'), ('CP001217', 'HPP12_0269')),
     (('AE000511', 'HP_1107'), ('CP001217', 'HPP12_1072')), (('AE000511', 'HP_0793'), ('CP001217', 'HPP12_0800')),
     (('AE000511', 'HP_1072'), ('CP001217', 'HPP12_0372')), (('AE000511', 'HP_0415'), ('CP001217', 'HPP12_1009')),
     (('AE000511', 'HP_1100'), ('CP001217', 'HPP12_1065')), (('AE000511', 'HP_0965'), ('CP001217', 'HPP12_0960')),
     (('AE000511', 'HP_0762'), ('CP001217', 'HPP12_0772'))}


    humanName = {}
    humanName['AE000511'] = '26695'
    humanName['CP001217'] = 'P12'
    humanName['CP001173'] = 'G27'
    humanName['AE001439'] = 'J99'

    intGenes = defaultdict(set)
    for p1, p2 in diffGenePairs:

        intGenes[p1[0]].add(p1[1])
        intGenes[p2[0]].add(p2[1])



    #highInterestGenes = {'HP_0965', 'HP_1141', 'HP_0099', 'HP_0064', 'HP_0818', 'HP_0312', 'HP_0108', 'HP_0342', 'HP_0235',                'HP_1107', 'HP_0104', 'HP_0860', 'HP_0656', 'HP_1533', 'HP_0963', 'HP_0755', 'HP_1184', 'HP_0142',                'HP_0762', 'HP_1455', 'HP_0043', 'HP_0783', 'HP_0744_2', 'HP_0599', 'HP_0759', 'HP_0743', 'HP_0091',                'HP_1277', 'HP_1046', 'HP_0897', 'HP_0519', 'HP_0048', 'HP_1349', 'HP_0289', 'HP_1213', 'HP_0558',                'HP_1105', 'HP_0286', 'HP_1229', 'HP_0430', 'HP_1282', 'HP_0270', 'HP_0415', 'HP_0036', 'HP_0868',                'HP_0793', 'HP_0566', 'HP_0661', 'HP_1072', 'HP_0568', 'HP_1100'}
    #interestGenes = {'HP_0580', 'HP_0208_2', 'HP_1533', 'HP_0342', 'HP_0793', 'HP_0566', 'HP_0783', 'HP_0906', 'HP_0661', 'HP_0099', 'HP_1141', 'HP_1486', 'HP_1259', 'HP_1216', 'HP_0731', 'HP_0460', 'HP_0104', 'HP_1517', 'HP_0547', 'HP_0430', 'HP_0568', 'HP_0770', 'HP_1184', 'HP_0272', 'HP_0599', 'HP_0860', 'HP_0048', 'HP_0505', 'HP_1265', 'HP_0922', 'HP_0963', 'HP_0898', 'HP_0788', 'HP_0091', 'HP_0431', 'HP_1105', 'HP_1499', 'HP_0289', 'HP_1471', 'HP_1277', 'HP_0286', 'HP_0428', 'HP_1072', 'HP_0270', 'HP_1177', 'HP_1349', 'HP_0066', 'HP_1229', 'HP_0047', 'HP_0312', 'HP_0142', 'HP_0818', 'HP_1078', 'HP_0694', 'HP_1361', 'HP_0298', 'HP_1044', 'HP_1513', 'HP_0995', 'HP_1282', 'HP_0012', 'HP_0868', 'HP_0854', 'HP_0755', 'HP_0434', 'HP_0759', 'HP_0338', 'HP_0456', 'HP_0762', 'HP_0108', 'HP_0681', 'HP_0743', 'HP_1157', 'HP_1566', 'HP_0030', 'HP_0760', 'HP_1046', 'HP_0884', 'HP_0465', 'HP_0238', 'HP_1522_0', 'HP_0656', 'HP_0927', 'HP_0519', 'HP_0421', 'HP_0807', 'HP_0235', 'HP_0897', 'HP_0747', 'HP_0965', 'HP_0517', 'HP_0744_2', 'HP_0415', 'HP_0059', 'HP_0262', 'HP_1455', 'HP_1107', 'HP_0064', 'HP_0043', 'HP_0036', 'HP_0608', 'HP_0964', 'HP_1213', 'HP_0558', 'HP_1100', 'HP_0462'}
    #print(len(interestGenes))



    for genome in allowedOrgs:
        genomeDB.loadGenome(genome)

    resCounts = defaultdict(lambda: defaultdict(list))

    for homID in homolDB.homologies:
        homGroup = homolDB.get_homology_cluster(homID)
        homCluster = [(x, homGroup[x]) for x in homGroup]

        orgResults = []

        for org in allowedOrgs:
            orgT = [x for x in homCluster if x[0] == org]

            orgResults.append(orgT)

        geneOfInterest = False

        for clusterElems in orgResults:

            for foundElem in clusterElems:

                if foundElem[0] in intGenes and foundElem[1] in intGenes[foundElem[0]]:
                    geneOfInterest = True


        #if len([x for x in orgResults if len(x) == 0]):
        #    continue

        orgT = [x[0] for x in orgResults if len(x) > 0]

        for idx, orgI in enumerate(orgT):
            seq = genomeDB.get_sequence( orgI[0], orgI[1] )

            aaCount = 0
            for aa in allAA:
                aaCount += seq.count(aa)

            relCounts = aaCount / len(seq)

            if geneOfInterest:
                resCounts[orgI[0]]['Interest'].append(relCounts)
            else:
                resCounts[orgI[0]]['Other'].append(relCounts)

    if False:

        for genome in allowedOrgs:

            resCounts = defaultdict(list)

            highInterestGenes = set()
            interestGenes = intGenes[genome]


            for gene in genomeDB.genomes[genome]:

                protSeq = genomeDB.get_sequence(genome, gene)

                aaCount = 0

                for aa in allAA:
                    aaCount += protSeq.count(aa)

                relCounts = aaCount / len(protSeq)

                if gene in highInterestGenes:
                    resCounts['Interest'].append(relCounts)
                    resCounts['High Interest'].append(relCounts)
                elif gene in interestGenes:
                    resCounts['Interest'].append(relCounts)
                else:
                    resCounts['Other'].append(relCounts)


    plt.style.use('seaborn-pastel')


    matplotlib.rcParams['font.size'] = 22
    matplotlib.rcParams['patch.linewidth'] =3


    def makeOrganismComparison(resCounts, genomeI, genomeJ, data='Interest'):

        fig, ax = plt.subplots(figsize=(9, 6))

        interestData = resCounts[genomeI][data]
        genomeData = resCounts[genomeJ][data]

        ksTestInt = ks_2samp(interestData, genomeData)

        histData = OrderedDict()

        histData[humanName.get(genomeI, genomeI) + " n=" + str(len(interestData)) + " "+data]= interestData
        histData[humanName.get(genomeJ, genomeJ) + " n=" + str(len(genomeData)) + " " + data]= genomeData

        histDataFlat = [histData[x] for x in histData]

        ax.hist(histDataFlat, 100, histtype='step', cumulative=1, normed=1, stacked=False)

        ax.legend([x for x in histData], loc='lower right')

        formatBase = "{0:.6G}" if ksTestInt.pvalue > 0.05 else "{0:.6e}"
        roundedPVal = formatBase.format(ksTestInt.pvalue)

        ax.set_xlim((0, 0.13))
        ax.set_title(
            "Distribution of " + ", ".join(allAA) + " in "+data.lower()+" proteins\n(KS-test pVal=" + roundedPVal + ")")
        plt.tight_layout()
        plt.savefig(fileLocation + "/plots/" + genomeI + "_" + genomeJ + "_comparison_"+data.lower()+".png")



    for genome in allowedOrgs:
        fig, ax = plt.subplots(figsize=(9,6))
        interestData = resCounts[genome]['Interest']
        genomeData = resCounts[genome]['Other']

        ksTestInt = ks_2samp(interestData, genomeData)

        histData = OrderedDict()

        histData[humanName.get(genome, genome) + " n=" + str(len(interestData)) + " Interest"] = interestData
        histData[humanName.get(genome, genome) + " n=" + str(len(genomeData)) + " Other"] = genomeData

        histDataFlat = [histData[x] for x in histData]

        ax.hist( histDataFlat, 100, histtype='step', cumulative=1, normed = 1, stacked = False )
        ax.legend([x for x in histData], loc='lower right')

        formatBase = "{0:.6G}" if ksTestInt.pvalue > 0.05 else "{0:.6e}"
        roundedPVal = formatBase.format(ksTestInt.pvalue)

        ax.set_title("Distribution of " + ", ".join(allAA) + " in " + humanName.get(genome, genome) + "\n(KS-test pVal=" + roundedPVal + ")")
        #plt.show()
        plt.tight_layout()
        plt.savefig(fileLocation + "/plots/" + genome + "_comparison.png")


    for i in range(0, len(allowedOrgs)):
        for j in range(i+1, len(allowedOrgs)):

            genomeI = allowedOrgs[i]
            genomeJ = allowedOrgs[j]

            makeOrganismComparison(resCounts, genomeI, genomeJ, 'Interest')
            makeOrganismComparison(resCounts, genomeI, genomeJ, 'Other')

            #plt.show()