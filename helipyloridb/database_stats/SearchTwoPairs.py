import sys,os

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.svm import LinearSVC

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation
from porestat.utils.DataFrame import DataFrame, DataRow

import pandas as pd
import seaborn as sns
import numpy as np

from sklearn import decomposition
from sklearn import datasets

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"
    fileLocation = "/home/mjoppich/ownCloud/data/hpyloriDB/"
    fileLocation = "/mnt/c/dev/data/haas/homdb/"


    homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp_split")
    genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)

    """
    for combid in homDB.combinations:
        elems = homDB.combinations[combid]

        homDB.homologies[combid] = elems

    homDB.finalize()
    homDB.save_to_file(fileLocation + "combed")
    """

    for orgname in homDB.get_all_organisms():
        genomDB.loadGenome(orgname)
    allorgs = list(homDB.get_all_organisms())

    allorgs.remove("6_Tx30a")


    allData = DataFrame()


    allData.addColumns(allorgs)
    homClusterIDs = []

    extra = ['AE001439', 'CP009259']
    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc and not x in extra and not x.startswith("6_")] # and not x.startswith("15")

    allorgs = extra + mc + nmc

    print("MC", len(mc), mc)
    print("NMC", len(nmc), nmc)


    for homID in homDB.homologies:

        clust = homDB.get_cluster(homID)


        val = homDB.get_homology_cluster(homID)
        maxlength = 0
        for org in val:

            geneid = val[org]
            seq = genomDB.get_sequence(org, geneid)

            if len(seq) > maxlength:
                maxlength = len(seq)


        if maxlength < 80:
            continue

        mcCount = sum([1 for x in mc if x in clust])
        nmcCount = sum([1 for x in nmc if x in clust])

        hasit15 = "15_N1-024A2" in [x for x in nmc if x in clust]
        hasit12 = "12_N4-050C1" in [x for x in nmc if x in clust]
        hasit1 = "1_N1-024A1" in [x for x in mc if x in clust]
        hasit14 = "14_1-20A_UB64" in [x for x in mc if x in clust]

        hasitJ99 = "AE001439" in clust
        hasitSS1 = "CP009259" in clust

        if mcCount == len(mc) and nmcCount == len(nmc):
            continue

        if not hasitJ99 or not hasitSS1:
            continue

        if not mcCount == len(mc):
            continue

        homClusterIDs.append(homID)


    print(allorgs)
    print(len(homClusterIDs))

    print(homClusterIDs)

    orgMatrix = []

    for org in allorgs:

        orgRes = []
        for homID in homClusterIDs:
            val = org in homDB.get_homology_cluster(homID)

            if val:
                orgRes.append(1)
            else:
                orgRes.append(0)

        orgMatrix.append(orgRes)


    import scipy
    import pylab
    import scipy.cluster.hierarchy as sch
    import seaborn as sns
    import pandas as pd

    # Generate random features and distance matrix.
    D = np.matrix(orgMatrix)

    Dt = D.transpose()

    df = pd.DataFrame(D, columns=homClusterIDs, index=allorgs)

    for name, values in df.iteritems():
        print('{name}: {value}'.format(name=name, value=values[0]))

    g = sns.clustermap(df, method="average",row_cluster=False )

    plt.show()
