import sys,os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase

if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp_split")
    genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)

    allorgs = homDB.get_all_organisms()

    extra = ['AE001439', 'CP009259']
    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc and not x in extra and not x.startswith("6_")] # and not x.startswith("15")

    print("MC", len(mc), mc)
    print("NMC", len(nmc), nmc)

    homlist = []

    for homID in homDB.homologies:

        clust = homDB.get_cluster(homID)

        mcCount = sum([1 for x in mc if x in clust])
        nmcCount = sum([1 for x in nmc if x in clust])

        hasit15 = "15_N1-024A2" in [x for x in nmc if x in clust]
        hasit12 = "12_N4-050C1" in [x for x in nmc if x in clust]
        hasit1 = "1_N1-024A1" in [x for x in mc if x in clust]
        hasit14 = "14_1-20A_UB64" in [x for x in mc if x in clust]

        hasitJ99 = "AE001439" in clust
        hasitSS1 = "CP009259" in clust

        if hasitJ99 and hasitSS1 and nmcCount == 0 and mcCount > 0:
            print(homID, mcCount, nmcCount, clust)
            homlist.append(homID)

        if not hasitSS1 and not hasitJ99 and mcCount == 0 and nmcCount > 0:
            print(homID, mcCount, nmcCount, clust)
            homlist.append(homID)

    print(homlist)


