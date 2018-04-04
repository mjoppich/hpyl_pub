import argparse
import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from analysis.homologybuilder import HomologyBuilder


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=str, help='input', required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='input', required=True)

    args = parser.parse_args()

    fileLocation = args.location

    builder = HomologyBuilder(basePath=fileLocation)
    homolDB = builder.analyse()
    homolDB.save_to_file(args.output.name)