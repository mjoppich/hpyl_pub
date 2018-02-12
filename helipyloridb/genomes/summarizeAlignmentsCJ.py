from analysis.homologybuilder import HomologyBuilder
from database.genomedb import GenomeDB

if __name__ == '__main__':
    fileLocation = '/mnt/c/ownCloud/data/cbdb/'


    initialise=False
    if initialise:

        genomDB = GenomeDB(fileLocation + "genomes/", fileFormat='gb', fileExtension='.gbff', loadAll=True)
        genomDB.writeBLASTfastas(fileLocation + "genomes/")

        exit()

    builder = HomologyBuilder(basePath=fileLocation, inputFormat='gb', inputExtension='.gbff')
    homolDB = builder.analyse()
    homolDB.save_to_file(fileLocation + "/cbj")