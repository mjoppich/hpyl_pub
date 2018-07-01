from collections import defaultdict

import os

from utils import fileLocation
from utils.DataFrame import DataFrame
from utils.GeneOntology import GeneOntology


class XRefDatabase:

    def __init__(self, fileName=fileLocation + "/hpp12_hp_xref"):

        self.df = DataFrame.parseFromFile(fileName)

        self.infos = {}
        self.add_infos = defaultdict(list)
        self.xrefs = {
            'Uniprot': 'GeneIdentity.UNIPROT',
        'GO': 'GeneIdentity.GO_ID',
        'Pfam': 'GeneIdentity.PFAM',
        'Interpro': 'GeneIdentity.INTERPRO',

        }

        self.go = GeneOntology("/mnt/c/ownCloud/data/" + "miRExplore/go/go.obo")

        for row in self.df:
            elemName = row['GeneIdentity.GENE_NAME']
            self.infos[elemName] = row


        if os.path.exists(fileName + "_add"):
            self.df_add = DataFrame.parseFromFile(fileName + "_add")
            for row in self.df_add:
                self.add_infos[row['XREF']].append(row['GOID'])


    def get_infos(self, elemName, default=None):

        if not elemName in self.infos:
            return default

        return self.infos[elemName]

    def make_infos(self, elemName):

        retElem = defaultdict(list)

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

        retElem['TMGO'] = []
        for x in self.xrefs:
            for vals in retElem[x]:
                if vals in self.add_infos:
                    for addI in self.add_infos[vals]:

                        if self.goNotIncluded(addI, retElem):
                            retElem['TMGO'].append(addI)


        retElem = self.addGOCats(retElem)

        retElem = self.restructureGO(retElem)

        return retElem

    def restructureGO(self, retElem):

        restructuredGO = {}
        restructuredGO['Uniprot'] = retElem['Uniprot']

        restructuredGO['Interpro'] = retElem['Interpro']
        restructuredGO['Pfam'] =retElem['Pfam']
        restructuredGO['GOTERMS'] = retElem['GOTERMS']

        restructuredGO['GO'] = defaultdict(list)
        restructuredGO['TMGO'] = defaultdict(list)

        for goID in retElem['GO']:
            exGO = self.go.getID(goID)
            for x in exGO.namespace:
                restructuredGO['GO'][x].append(goID)

        for goID in retElem['TMGO']:
            exGO = self.go.getID(goID)
            for x in exGO.namespace:
                restructuredGO['TMGO'][x].append(goID)

        return restructuredGO

    def addGOCats(self, retElems):

        retElems['GOTERMS'] = {}

        for goID in retElems['GO'] + retElems['TMGO']:
            exGO = self.go.getID(goID)

            retElems['GOTERMS'][goID] = exGO.name

        return retElems

    def goNotIncluded(self, goID, retElems):

        if goID in retElems['GO']:
            return False

        if goID in retElems['TMGO']:
            return False

        for exGoID in retElems['GO'] + retElems['TMGO']:
            exGO = self.go.getID(exGoID)

            if goID in [x.termid for x in exGO.getAllChildren()]:
                return False

        return True


if __name__ == '__main__':

    db = XRefDatabase()

    row = db.get_infos('jhp_1366')
    print(row.to_pairs())

    print(db.make_infos('jhp_1366'))