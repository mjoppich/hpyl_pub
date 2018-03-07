from enum import Enum

class GeneIdentity(Enum):
    GENE_SYMBOL = 'symbol'
    ALIAS = "alias_symbol"
    PREV_NAME = "prev_name"
    ORDERED_LOCUS = 'ordered_locus'
    ORF_NAME = 'orf_name'

    ENTREZ = "entrez_id"
    ENSEMBL = "ensembl_gene_id"
    VEGA = "vega_id"
    VIRUS_HOST = 'virus_host'
    GENE_NUMBER='P_GI'


    UNIPROT = "uniprot_ids"
    UNIPROT_NAME= 'entry name'

    HGNC_ID = "hgnc_id"
    MGD_ID = "mgd_id"

    ORGANISM="organism"
    TAXID = "taxid"
    TAX_LINEAGE = 'lineage'

    GENE_FAMILY = "gene_family"
    STATUS = "status"

    PROTEOMES="proteome"
    PROTEIN_NAMES="protein names"

    ACCESSION='accession'
    SEQ_TYPE='seq_type'
    SEQUENCE='sequence'
    SEQUENCE_LENGTH='sequence_length'
    GENE_NAME='GENENAME'
    GO='GO'
    GO_ID='GO_ID'
    PFAM='families'
    INTERPRO='interpro'

    def describe(self):
        # self is the member here
        return self.name, self.value


class GeneEntityException(Exception):

    def __init__(self, msg):
        super(GeneEntityException, self).__init__()
        self.msg = msg

class GeneEntity:

    def __init__(self, props):

        self.infos = props

    def get(self, entity, default = None):

        if type(entity) == GeneIdentity:
            entity = [entity]

        fetchEntites = []

        for x in entity:
            if x in self.infos:
                fetchEntites.append(x)

        if len(fetchEntites) == 0:
            raise GeneEntityException("no known entity included")

        dReturn = {}

        for x in fetchEntites:
            dReturn[x] = self.infos[x]

        return dReturn