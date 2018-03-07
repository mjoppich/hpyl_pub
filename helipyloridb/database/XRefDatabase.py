from collections import defaultdict
from porestat.utils.DataFrame import DataFrame

from utils import fileLocation


class XRefDatabase:

    def __init__(self, fileName=fileLocation + "/hpp12_hp_xref"):

        self.df = DataFrame.parseFromFile(fileName)

        self.infos = {}
        self.xrefs = {
            'Uniprot': 'GeneIdentity.UNIPROT',
        'GO': 'GeneIdentity.GO_ID',
        'Pfam': 'GeneIdentity.PFAM',
        'Interpro': 'GeneIdentity.INTERPRO',

        }

        for row in self.df:
            elemName = row['GeneIdentity.GENE_NAME']
            self.infos[elemName] = row


    def get_infos(self, elemName, default=None):

        if not elemName in self.infos:
            return default

        return self.infos[elemName]

    def make_infos(self, elemName):

        retElem = {}

        for x in self.xrefs:
            retElem[x] = []

        elemInfo = self.get_infos(elemName)

        if elemInfo == None:
            return retElem

        for x in self.xrefs:
            info = elemInfo[self.xrefs[x]]

            if info == None or info == 'None':
                continue

            for val in info.split(';'):
                vals = val.strip()

                if len(vals) == 0:
                    continue
                retElem[x].append(vals)

        return retElem



if __name__ == '__main__':

    db = XRefDatabase()

    row = db.get_infos('jhp_1366')
    print(row.to_pairs())

    print(db.make_infos('jhp_1366'))