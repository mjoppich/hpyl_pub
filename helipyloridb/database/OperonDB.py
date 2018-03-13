from collections import defaultdict

from openpyxl import load_workbook


class OperonDB:

    def __init__(self):

        self.operons = defaultdict(list)
        self.operon_props = defaultdict(lambda: dict())


    def find_gene(self, geneID, defaultval=list()):

        retOperonIDs = []

        for x in self.operons:

            if geneID in self.operons[x]:
                retOperonIDs.append(x)

        return defaultval if len(retOperonIDs) == 0 else retOperonIDs

    def get_operon_infos(self, operonID):

        if not operonID in self.operons:
            return None


        operonBase = {
            'OPERONID': operonID,
            'GENES': self.operons[operonID]
        }

        operonProps = self.operon_props[operonID]

        for x in operonProps:
            operonBase[x] = operonProps[x]

        return operonBase


    @classmethod
    def from_cs_operons(cls, filelocation="/mnt/c/ownCloud/data/hpyloriDB/sharma/operons.xlsx"):

        opDB = OperonDB()

        wb = load_workbook(filelocation)
        ws = wb.active

        curOperonID = 0

        for row in ws:
            if row[0].row < 7:
                continue

            doorIdx = str(row[0].value)

            if doorIdx != None and doorIdx != "None":
                doorIdx = doorIdx.split(", ")
            else:
                doorIdx = []

            strand = row[1].value
            evaluationResult = row[2].value
            comment = row[3].value

            genes = []
            for i in range(4, len(row)):
                rowVal = row[i].value

                if rowVal == None:
                    continue

                genes.append(rowVal)

            modGenes = [x[0:2]+'_' + x[2:] for x in genes]

            curOpID = 'OPERON' + str(curOperonID)

            opDB.operon_props[curOpID]['strand'] = strand
            opDB.operon_props[curOpID]['evaluation'] = evaluationResult
            opDB.operon_props[curOpID]['comment'] = comment if comment != None else ""
            opDB.operon_props[curOpID]['strand'] = strand
            opDB.operon_props[curOpID]['DOOR'] = doorIdx

            opDB.operons[curOpID] = modGenes

            curOperonID += 1

            #print(curOpID, doorIdx, strand, evaluationResult, comment, modGenes)
        return opDB


if __name__ == '__main__':

    OperonDB.from_cs_operons()
