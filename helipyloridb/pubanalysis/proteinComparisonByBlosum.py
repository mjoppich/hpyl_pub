import sys
import os
from collections import Counter, defaultdict

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

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

if __name__ == '__main__':

    #compareAA = (compareAA[0], [x for x in aaUsage['CP001217'] if x != '*' and not x in compareAA[0]])


    def calculateDifferences(homolDB, orgI, orgJ,):

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

            if '*' in seqI or '*' in seqJ:
                continue

            if 'X' in seqI or 'X' in seqJ:
                continue

            if 'U' in seqI or 'U' in seqJ:
                continue

            print("HOMID", orgI, orgJ, homID)


            aaLength += (len(seqI)+len(seqJ))/2.0
            alignments = pairwise2.align.globalds(seqI, seqJ, matrix, -50, -.2)

            alignment = alignments[0]

            alignI = alignment[0]
            alignJ = alignment[1]

            for i in range(0, len(alignI)):

                aaI = alignI[i]
                aaJ = alignJ[i]

                if aaI == '-':
                    orgSubMatrix[(aaI, aaJ)] += 1
                elif aaJ == '-':
                    orgSubMatrix[(aaI, aaJ)] += 1
                else:

                    if aaI < aaJ:
                        orgSubMatrix[(aaI, aaJ)] += 1
                    else:
                        orgSubMatrix[(aaJ, aaI)] += 1


        for x in orgSubMatrix:
            orgSubMatRel[x] = orgSubMatrix[x] / aaLength

        return orgSubMatRel


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

    print("Starting Differences")
    (orgI, orgJ) = ('NC003912', 'NC002163')
    subMatrix[(orgI, orgJ)] = calculateDifferences(cbHomolDB, orgI, orgJ)


    print("Starting HPP")
    (orgI, orgJ) = ('CP001217', 'AE000511')
    subMatrix[(orgI, orgJ)] = calculateDifferences(hpHomolDB, orgI, orgJ)

    print(subMatrix)