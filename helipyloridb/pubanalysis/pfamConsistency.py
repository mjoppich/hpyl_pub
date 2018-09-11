import argparse
import os
from collections import Counter

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase

import matplotlib.pyplot as plt

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=argparse.FileType('r'), help='input', required=True)


    args = parser.parse_args()

    homDB = HomologyDatabase.loadFromFile(args.location.name)

    genomDB = GenomeDB(os.path.dirname(args.location.name) + "/genomes", loadAll=False)
    analyse = HomDBAnalyser(homDB, genomDB)


    sizeCounter = Counter()
    seqsPerCluster = []
    errors = []

    print("Performing calcs for cluster count", len(homDB.homologies))


    for idx, homID in enumerate(homDB.homologies):

