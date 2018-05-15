import argparse
import sys, os

from database.homologydb import HomologyDatabase

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from analysis.homologybuilder import HomologyBuilder


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=argparse.FileType('r'), help='input', required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='input', required=True)

    args = parser.parse_args()

    homDB = HomologyDatabase.loadFromFile(args.location.name)

    for combID in homDB.combinations:

        cluster = homDB.combinations[combID]
        cluster = list(cluster)

        for i in range(0, len(cluster)):

            for j in range(i+1, len(cluster)):

                homDB.addHomologyRelation( cluster[i], cluster[j], None )


    args.output.close()

    homDB.finalize()
    homDB.save_to_file(args.output.name)