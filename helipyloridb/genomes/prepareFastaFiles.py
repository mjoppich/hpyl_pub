import glob
from Bio import SeqIO
import argparse

import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from database.genomedb import GenomeDB



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=str, help='input', required=True)
    args = parser.parse_args()

    fileLocation = args.location

    for file in glob.glob(fileLocation + '/*.gb'):

        print(file)

        genomeDB = GenomeDB(fileLocation, loadAll=False)
        genomeDB.loadGenome(file, True)

        genomeDB.writeBLASTfastas(fileLocation)

