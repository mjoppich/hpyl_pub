import os
from collections import defaultdict

from database.DiamondResult import DiamondResult
from database.genomedb import GenomeDB


class GeneDuplicationDB:

    def __init__(self):

        #org -> genesets
        self.duplicatedGenes = defaultdict( list )

    def __str__(self):

        seq = ""
        for org in self.duplicatedGenes:
            for dupSet in self.duplicatedGenes[org]:
                seq += "\t".join([str(org)] + [str(x) for x in dupSet]) + "\n"

        return seq

    def print(self):

        for org in self.duplicatedGenes:
            for dupSet in self.duplicatedGenes[org]:
                print("\t".join([str(org)] + [str(x) for x in dupSet]))


    def load_organism(self, fp, orgGenomeDB=None):

        with open(fp, 'r') as infile:

            genomeID = str(os.path.basename(fp).split(".")[0])

            if orgGenomeDB == None:
                orgGenomeDB = GenomeDB(os.path.dirname(fp) + "../genomes/"+genomeID + ".fa")

            for line in infile:
                ret = DiamondResult.from_line(line, genomeID, genomeID)

                if ret.identity < 0.95:
                    continue

                if ret.subject.seqid == ret.query.seqid:
                    continue

                subjSeq = orgGenomeDB.get_sequence(genomeID, ret.subject.seqid)
                querySeq = orgGenomeDB.get_sequence(genomeID, ret.query.seqid)

                partialSQ = (len(subjSeq)/len(querySeq))
                partialQS = (len(querySeq) / len(subjSeq))

                partialSQok = 0.95 < partialSQ and partialSQ < 1.05
                partialQSok = 0.95 < partialQS and partialQS < 1.05

                if not partialQSok and not partialSQok:
                    continue

                self.add_gene_duplication(genomeID, ret.subject.seqid, ret.query.seqid)




    def add_gene_duplication(self, organism, gene1, gene2):

        orgGenes = self.duplicatedGenes[organism]

        for i in range(0, len(orgGenes)):

            dupset = orgGenes[i]

            if gene1 in dupset or gene2 in dupset:
                dupset.append(gene1)
                dupset.append(gene2)

                orgGenes[i] = sorted(list(set(dupset)))

                return


        newdupset = list()
        newdupset.append(gene1)
        newdupset.append(gene2)

        orgGenes.append( sorted(newdupset) )

        return

    def get_gene_duplication(self, organism, gene):

        orgGenes = self.duplicatedGenes.get(organism, None)

        if orgGenes == None:
            return None

        for i in range(0, len(orgGenes)):

            dupset = orgGenes[i]

            if gene in dupset:
                return dupset

        return None

    def has_gene_duplication(self, organism, gene):

        return self.get_gene_duplication(organism, gene) != None
