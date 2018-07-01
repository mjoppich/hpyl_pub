import sys,os

from database.homDBAnalyser import HomDBAnalyser

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
from porestat.utils.DataFrame import DataFrame, DataRow

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def intOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def distance(r1, r2):
    # sort the two ranges such that the range with smaller first element
    # is assigned to x and the bigger one is assigned to y
    x, y = sorted((r1, r2))

    if x[0] <= x[1] < y[0] and all(y[0] <= y[1] for y in (r1, r2)):
        return y[0] - x[1]
    return 0

if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp_split")
    genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)

    allorgs = homDB.get_all_organisms()

    for org in allorgs:
        genomDB.loadGenome(org)

    extra = ['AE001439', 'CP009259']
    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc and not x in extra and not x.startswith("6_")] # and not x.startswith("15")

    print("MC", len(mc), mc)
    print("NMC", len(nmc), nmc)

    analyse = HomDBAnalyser(homDB, genomDB)


    for homID in homDB.homologies:

        homCluster = homDB.homologies[homID]


        mcCount = 0
        nmcCount = 0

        org2length = {}

        for (org, geneid) in homCluster:
            if org in mc:
                mcCount += 1

            if org in nmc:
                nmcCount += 1

            seq = genomDB.get_sequence(org, geneid)

            org2length[org] = seq


        mcLengths = set()
        nmcLenghts = set()

        mcSeqs = set()
        nmcSeqs = set()

        for org in org2length:

            if org in mc:
                seq = org2length[org]
                mcLengths.add(len(seq))
                mcSeqs.add(seq)

            if org in nmc:
                seq = org2length[org]
                nmcLenghts.add(len(seq))
                nmcSeqs.add(seq)

        if len(mcLengths) == 0 or len(nmcLenghts) == 0:
            continue


        mcElemPair = (min(mcLengths), max(mcLengths))
        nmcElemPair = (min(nmcLenghts), max(nmcLenghts))

        maxDistance = distance(nmcElemPair, mcElemPair)
        maxOverlap = intOverlap(mcElemPair, nmcElemPair)

        if maxDistance > 10 or nmcCount > 1.5*len(nmc):
            print(homID, mcCount, nmcCount, maxDistance, mcElemPair, nmcElemPair)
            #print(mcSeqs)
            #print(nmcSeqs)


