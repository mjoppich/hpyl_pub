from analysis.homologybuilder import HomologyBuilder
from utils import fileLocation


if __name__ == '__main__':

    builder = HomologyBuilder(basePath=fileLocation)
    homolDB = builder.analyse()
    homolDB.save_to_file(fileLocation + "/hpp12_hp")