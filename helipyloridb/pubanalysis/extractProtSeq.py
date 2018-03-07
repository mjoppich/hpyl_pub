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


    genomeDB = GenomeDB(fileLocation + "/genomes/")

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

    allGenomes = set([x[0][0] for x in diffGenePairs]).union(set([x[1][0] for x in diffGenePairs]))

    for genome in allGenomes:
        genomeDB.loadGenome(genome)


    fastaLines = ([], [])

    for x1, x2 in diffGenePairs:
        seq1 = genomeDB.get_sequence(x1[0], x1[1])
        seq2 = genomeDB.get_sequence(x2[0], x2[1])

        if seq1 == None:
            print(x1)

        if seq2 == None:
            print(x2)

        fastaLines[0].append(">" + x1[1] + " " + x1[0] + "\n" + seq1)
        fastaLines[1].append(">" + x2[1] + " " + x2[0] + "\n" + seq2)


    print("\n".join(fastaLines[0]))

    print()
    print()
    print()

    print("\n".join(fastaLines[1]))