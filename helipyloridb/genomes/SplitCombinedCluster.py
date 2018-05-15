import argparse
import sys, os

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from analysis.homologybuilder import HomologyBuilder


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=argparse.FileType('r'), help='input', required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='input', required=True)

    args = parser.parse_args()

    homDB = HomologyDatabase.loadFromFile(args.location.name)

    genomDB = GenomeDB(os.path.dirname(args.location.name) + "/genomes", loadAll=False)

    allorgs = homDB.get_all_organisms()

    for org in allorgs:
        genomDB.loadGenome(org)

    analyse = HomDBAnalyser(homDB, genomDB)


    for homID in homDB.homologies:

        if not homID == 'HOMID139':
            continue

        cluster = homDB.homologies[homID]

        if len(cluster) > 30:

            print(homID, len(cluster))



            aligned = analyse.cluster_align_clustalw(homID)
            for rec in sorted(aligned._records, key=lambda x: x.id):

                print(rec.seq, rec.id)


            calculator = DistanceCalculator('identity')
            constructor = DistanceTreeConstructor(calculator, 'upgma')
            tree = constructor.build_tree(aligned)
            dm = calculator.get_distance(aligned)
            tree.ladderize()  # Flip branches so deeper clades are displayed at top

            Phylo.draw_ascii(tree)
            Phylo.draw(tree)

            print(dm)





    #homDB.save_to_file(args.output.name)