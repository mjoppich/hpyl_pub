import glob
from collections import defaultdict

from Bio import SeqIO


class GenomeDB:

    def __init__(self, location, loadAll = False):

        self.genomes = defaultdict(lambda: dict())

        if loadAll:
            for file in glob.glob(location+'/*.gb'):
                self.loadGenome(file)

    def loadGenome(self, file):

        gbParser = SeqIO.parse(file, "embl")

        for gb_record in gbParser:

            genomeID = gb_record.name

            for feature in gb_record.features:

                if not feature.type.upper() == 'CDS':
                    continue

                locTag = feature.qualifiers.get('locus_tag', [None])[0]
                proID = feature.qualifiers.get('protein_id', [None])[0]
                translation = feature.qualifiers.get('translation', [None])[0]

                productID = locTag if locTag != None else proID

                if productID == None:
                    print("CDS without id:")
                    print(feature)
                    continue

                if translation == None:
                    continue

                self.genomes[genomeID][productID] = translation

        print("Loaded Genome: " + file)

    def get_sequence(self, genome, productID):
        return self.genomes.get(genome, {}).get(productID, None)