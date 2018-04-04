
from io import StringIO

from porestat.utils.DataFrame import DataFrame, DataRow

from xrefs.GeneIdentity import GeneIdentity
from xrefs.Store import RESTStore, StoreException, RequestMethod
import json

class UniprotStore(RESTStore):

    def __init__(self):

        super(UniprotStore, self).__init__( "Uniprot", "http://www.uniprot.org/uploadlists/" )

        self.columns = {
            GeneIdentity.UNIPROT: 'id',
            GeneIdentity.UNIPROT_NAME: 'entry name',
            GeneIdentity.ALIAS: 'genes',
            GeneIdentity.GENE_SYMBOL: 'genes(PREFERRED)',
            GeneIdentity.ORDERED_LOCUS: 'genes(OLN)',
            GeneIdentity.ORF_NAME: 'genes(ORF)',
            GeneIdentity.ORGANISM: 'organism',
            GeneIdentity.TAXID: 'organism-id',
            GeneIdentity.PROTEIN_NAMES: 'protein names',
            GeneIdentity.PROTEOMES: 'proteome',
            GeneIdentity.TAX_LINEAGE: 'lineage(all)',
            GeneIdentity.VIRUS_HOST: 'virus host',
            GeneIdentity.GO: 'go',
            GeneIdentity.GO_ID: 'go-id',
            GeneIdentity.SEQUENCE_LENGTH: 'length',
            GeneIdentity.PFAM: 'database(Pfam)',
            GeneIdentity.INTERPRO: 'database(Interpro)'

        }

        self.fromEntities = {
            GeneIdentity.UNIPROT: 'ID',
            GeneIdentity.ENSEMBL: 'EMBL_ID',
            GeneIdentity.GENE_NUMBER: 'P_GI',
            GeneIdentity.ENTREZ: 'P_ENTREZGENEID',
            GeneIdentity.GENE_NAME: 'GENENAME'
        }

        self.provides = [
            x for x in self.columns
        ]

        self.accepts = [
            x for x in self.fromEntities
        ]


    def _check_accepted(self, tocheck):

        if not type(tocheck) == list:
            tocheck = [tocheck]

        for x in tocheck:
            if not x in self.accepts:
                return False

        return True

    def _must_accept(self, tocheck):

        if not self._check_accepted(tocheck):
            raise StoreException("Not all entities are accepted by this store! \n Accepted %s \n Requested %s" % (self.accepts, tocheck))

        return True

    def _check_provided(self, tocheck):

        if not type(tocheck) == list:
            tocheck = [tocheck]

        for x in tocheck:
            if not x in self.provides:
                return False

        return True

    def _must_provide(self, tocheck):

        if not self._check_provided(tocheck):
            raise StoreException("Not all entities are provided by this store! \n Provided %s \n Requested %s" % (self.provides, tocheck))

        return True

    def _make_params(self, queryType, elements, toType):

        reqParams = {}

        reqParams['from'] = self.fromEntities[queryType]
        reqParams['to'] = 'ACC'

        reqParams['format'] = 'tab'
        reqParams['columns'] = ",".join([self.columns[x] for x in toType])

        reqParams['query'] = " ".join([str(x) for x in elements])

        return reqParams



    def fetch(self, fromEntity, elements, toEntities = [GeneIdentity.UNIPROT, GeneIdentity.GENE_SYMBOL, GeneIdentity.ORDERED_LOCUS], error_on_empty_result=True):

        self._must_accept(fromEntity)
        self._must_provide(toEntities)

        elements = sorted(elements)

        reqParams = self._make_params(fromEntity, elements, toEntities)

        for x in reqParams:

            lenReqParams = len(reqParams[x])

            if lenReqParams < 100:
                print(str(x) + " " + str(reqParams[x]))
            else:
                print(str(x) + " " + str(lenReqParams) + " elements")

        resp = self._request( RequestMethod.POST, "", reqParams)

        if (resp.text == None):
            print(json.dumps(reqParams))
            raise StoreException("Could not retrieve elements")

        if len(resp.text) == 0 and error_on_empty_result:
            raise StoreException("Empty result")

        convData = DataFrame()
        dfCols = toEntities + [fromEntity]
        convData.addColumns(dfCols)


        def addLineToReturn(lineData):

            modLine = {}
            for c,x in zip(dfCols, lineData):
                if x == '':
                    modLine[c] = None
                else:
                    modLine[c] = x

            convData.addRow( DataRow.fromDict(modLine) )


        bFirstLine = True
        for line in resp.text.split('\n'):

            if bFirstLine:
                bFirstLine = False
                continue

            if len(line) == 0:
                continue

            aline = line.split('\t')
            if len(aline) == 0:
                continue

            aline = aline[:-1]

            if ',' in aline[-1]:
                elems = aline[-1].split(',')
                elemCount = len(elems)

                for i in range (0, elemCount):

                    modLine = []

                    for elem in aline[:-1]:

                        aelem = elem.split(' ')
                        if len(aelem) != elemCount:
                            modLine.append(elem)
                        else:
                            modLine.append(aelem[i])

                    modLine.append(elems[i])

                    addLineToReturn(modLine)

            else:
                addLineToReturn(aline)

        return convData
