import sys
import os
from collections import Counter, defaultdict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils.utils import fileLocation

if __name__ == '__main__':

    #compareAA = (compareAA[0], [x for x in aaUsage['CP001217'] if x != '*' and not x in compareAA[0]])

    interestGenes = {'HP_0235', 'HP_1105', 'HP_1046', 'HP_0104', 'HP_0762', 'HP_0599', 'HP_0963', 'HP_0048', 'HP_1533', 'HP_0743', 'HP_0868', 'HP_0965', 'HP_0342', 'HP_0312', 'HP_0289', 'HP_0415', 'HP_0755', 'HP_0036', 'HP_1100', 'HP_1213', 'HP_0568', 'HP_0793', 'HP_1277', 'HP_1455', 'HP_1072', 'HP_1141', 'HP_0744_2', 'HP_0566', 'HP_0091', 'HP_1229', 'HP_0783', 'HP_0430', 'HP_1184', 'HP_0897', 'HP_0759', 'HP_1349', 'HP_1282', 'HP_0064', 'HP_0519', 'HP_0286', 'HP_0108', 'HP_0043', 'HP_0860', 'HP_0656', 'HP_1107', 'HP_0818', 'HP_0270', 'HP_0661', 'HP_0142', 'HP_0099', 'HP_0558'}
    #interestGenes = {'HP_0235', 'HP_0505', 'HP_0762', 'HP_0599', 'HP_0788', 'HP_1517', 'HP_1533', 'HP_0434', 'HP_0048', 'HP_0868', 'HP_0965', 'HP_0747', 'HP_0608', 'HP_0030', 'HP_0066', 'HP_1455', 'HP_1513', 'HP_0694', 'HP_1229', 'HP_0760', 'HP_0430', 'HP_0897', 'HP_0286', 'HP_0043', 'HP_1107', 'HP_1044', 'HP_0558', 'HP_1105', 'HP_0104', 'HP_0743', 'HP_1216', 'HP_1471', 'HP_1265', 'HP_0342', 'HP_0770', 'HP_0415', 'HP_0927', 'HP_0431', 'HP_1100', 'HP_1213', 'HP_0047', 'HP_1141', 'HP_0854', 'HP_0681', 'HP_0091', 'HP_0580', 'HP_0460', 'HP_1282', 'HP_0064', 'HP_1486', 'HP_0519', 'HP_1361', 'HP_1259', 'HP_1046', 'HP_0922', 'HP_1499', 'HP_0964', 'HP_0884', 'HP_0312', 'HP_0289', 'HP_0755', 'HP_0036', 'HP_1177', 'HP_0793', 'HP_0898', 'HP_1277', 'HP_1072', 'HP_0272', 'HP_0744_2', 'HP_0731', 'HP_0783', 'HP_1184', 'HP_0465', 'HP_1349', 'HP_0807', 'HP_0108', 'HP_0656', 'HP_0421', 'HP_0270', 'HP_0208_2', 'HP_0298', 'HP_1157', 'HP_0462', 'HP_0238', 'HP_1522_0', 'HP_0963', 'HP_0547', 'HP_0338', 'HP_0568', 'HP_1566', 'HP_1078', 'HP_0566', 'HP_0456', 'HP_0012', 'HP_0995', 'HP_0517', 'HP_0759', 'HP_0059', 'HP_0262', 'HP_0906', 'HP_0860', 'HP_0818', 'HP_0661', 'HP_0142', 'HP_0428', 'HP_0099'}


    def calculateDifferences(homolDB, orgI, orgJ,):

        orgSubMatrix = Counter()
        orgSubMatRel = Counter()
        orgSubMatrixDir = Counter()

        orgAACounts = defaultdict(lambda: Counter())

        interestCounter = 0

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

            if '*' in seqI or '*' in seqJ:
                continue

            if 'X' in seqI or 'X' in seqJ:
                continue

            if 'U' in seqI or 'U' in seqJ:
                continue

            #print("HOMID", orgI, orgJ, homID)


            aaLength += (len(seqI)+len(seqJ))/2.0
            alignments = pairwise2.align.globalds(seqI, seqJ, matrix, -50, -.2)

            alignment = alignments[0]

            alignI = alignment[0]
            alignJ = alignment[1]

            printAlignment=False

            for i in range(0, len(alignI)):

                aaI = alignI[i]
                aaJ = alignJ[i]

                orgSubMatrixDir[(aaI, aaJ)] += 1

                if set([aaI, aaJ]) == set(['W', 'L']):
                    if orgJT[1] in interestGenes:
                        interestCounter += 1
                        print(orgIT, orgJT, interestCounter)
                        #printAlignment = True

                    #printAlignment = False


                if aaI == '-':
                    orgSubMatrix[(aaI, aaJ)] += 1
                    orgAACounts[orgJT[0]][aaJ] += 1

                elif aaJ == '-':
                    orgSubMatrix[(aaI, aaJ)] += 1
                    orgAACounts[orgIT[0]][aaI] += 1

                else:

                    if aaI < aaJ:
                        pass#orgSubMatrix[(aaI, aaJ)] += 1
                    else:
                        pass
                    orgSubMatrix[(aaJ, aaI)] += 1

                    orgAACounts[orgIT[0]][aaI] += 1
                    orgAACounts[orgJT[0]][aaJ] += 1

            if printAlignment:
                print(orgJT[1])
                #print("WL Changes", orgIT, orgJT, orgSubMatrix[('L', 'W')])
                #print(alignment[0])
                #print(alignment[1])

        for x in orgSubMatrix:
            orgSubMatRel[x] = orgSubMatrix[x] / aaLength


        print("W content")
        for org in orgAACounts:
            print(org, "W", orgAACounts[org]['W'])

        return orgSubMatrix, orgSubMatRel, aaLength, orgSubMatrixDir


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
    (orgI, orgJ) = ('CP001217', 'AE000511')

    absMat, relMat, aaCount, dirMat = calculateDifferences(hpHomolDB, orgI, orgJ)

    print((orgI, orgJ))
    print(absMat)
    print(relMat)
    print(aaCount)
    print(dirMat)


    print("Starting Differences")
    (orgI, orgJ) = ('NC003912', 'NC002163')

    absMat, relMat, aaCount, dirMat = calculateDifferences(cbHomolDB, orgI, orgJ)

    print((orgI, orgJ))
    print(absMat)
    print(relMat)
    print(aaCount)
    print(dirMat)


