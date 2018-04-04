import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from analysis.homologybuilder import HomologyBuilder
from utils import fileLocation


if __name__ == '__main__':

    builder = HomologyBuilder(basePath=fileLocation)
    homolDB = builder.analyse()
    homolDB.save_to_file(fileLocation + "/hpdb_full")