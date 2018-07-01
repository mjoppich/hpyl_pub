import sys
import os
from collections import Counter

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils.utils import fileLocation

if __name__ == '__main__':


    homolDB = HomologyDatabase.loadFromFile(fileLocation+ "/hpp12_hp")
    genomeDB = GenomeDB(fileLocation + "/genomes/")

    allowedOrgs = ['CP001217', 'AE000511']

    compareAA = (['W'], ['F', 'G', 'A'])
    compareAA = (['W'], ['H', 'F', 'Y', 'P', 'K'])
    #compareAA = (['W', 'M'], ['H', 'F', 'Y', 'P', 'K'])
    #compareAA = (['W', 'M'], ['F', 'G', 'A'])



    for org in homolDB.get_all_organisms():
        if org in allowedOrgs:
            genomeDB.loadGenome(org)


    for org in allowedOrgs:

        orgAA = Counter()

        orgProts = genomeDB.genomes[org]
        totalLength = 0

        for protID in orgProts:

            AAseq = orgProts[protID]

            totalLength += len(AAseq)

            AAcounter = Counter(AAseq)

            orgAA += AAcounter

        print(org)
        print(orgAA)

        orgAArel = {}
        for aa in orgAA:
            orgAArel[aa] = orgAA[aa] / totalLength

        print("aaUsage[\'"+  org +"\'] = " + str(orgAArel))