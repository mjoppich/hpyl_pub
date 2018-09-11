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


    for idx, homID in enumerate(['HOMID95']):

        cluster = homDB.get_cluster(homID)
        sizeCounter[len(cluster)] += 1

        print("Organism count", len(cluster))
        print("Cluster Size", sum([len(x) for x in cluster]))

        mycluster = {}

        allowedOrgs = ['AE001439', 'AE000511', 'CP001217']
        for org in allowedOrgs:
            if org in cluster:
                mycluster[org] = cluster[org]

        print(mycluster)
        for org in mycluster:
            print(org, mycluster[org])


        aligned = analyse.cluster_align_clustalo(homID, allowedOrganisms=allowedOrgs)

        calculator = DistanceCalculator('blosum80')

        dm = calculator.get_distance(aligned)

        print(dm)

        pwError = sum([sum(x) for x in dm.matrix])
        errors.append(pwError)

        for elem in aligned:
            print(str(homID).rjust(5), elem.id.rjust(30), elem.seq)