from database.homologydb import HomologyDatabase
from utils import fileLocation
from xrefs.GeneIdentity import GeneIdentity
from xrefs.uniprotStore import UniprotStore

hpHomolDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp12_hp")

allIDS = []

for org in hpHomolDB.get_all_organisms():

    allOrgElems = hpHomolDB.get_organism_elements(org)

    for x in allOrgElems:
        allIDS.append(x)

up = UniprotStore()
allConvertedIDs = up.fetch(GeneIdentity.GENE_NAME, allIDS, toEntities=[GeneIdentity.UNIPROT, GeneIdentity.GO_ID, GeneIdentity.PFAM, GeneIdentity.INTERPRO, GeneIdentity.GENE_SYMBOL])

allConvertedIDs.export(fileLocation + "/hpp12_hp_xref")