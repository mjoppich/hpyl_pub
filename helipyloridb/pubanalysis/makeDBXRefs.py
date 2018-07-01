from database.homologydb import HomologyDatabase
from utils.utils import fileLocation
from xrefs.GeneIdentity import GeneIdentity
from xrefs.uniprotStore import UniprotStore

hpHomolDB = HomologyDatabase.loadFromFile(fileLocation + "/hpdb_full")

allIDS = []

for org in hpHomolDB.get_all_organisms():

    allOrgElems = hpHomolDB.get_organism_elements(org)

    for x in allOrgElems:
        allIDS.append(x)

print("Fetching information for", len(allIDS), "gene ids")


finalDF = None
chunksize = 1000

for i in range(0, len(allIDS), chunksize):

    imax = min([i+chunksize, len(allIDS)])

    chunkElems = allIDS[i:imax]

    print(i, imax, len(chunkElems))

    up = UniprotStore()
    allConvertedIDs = up.fetch(GeneIdentity.GENE_NAME, chunkElems, toEntities=[GeneIdentity.UNIPROT, GeneIdentity.GO_ID, GeneIdentity.PFAM, GeneIdentity.INTERPRO, GeneIdentity.GENE_SYMBOL], error_on_empty_result=False)

    if finalDF == None:
        finalDF = allConvertedIDs
    else:
        finalDF.merge(allConvertedIDs)

finalDF.export(fileLocation + "/hpdb_full_xref")