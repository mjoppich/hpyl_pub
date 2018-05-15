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

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpdb_full_new")
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


    allData = DataFrame()


    allData.addColumns(allorgs)
    homClusterIDs = []


    for homid in homDB.homologies:

        val = homDB.get_homology_cluster(homid)

        homClusterIDs.append(homid)

    print(allorgs)
    print(len(homClusterIDs))

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

    # Generate random features and distance matrix.
    D = np.matrix(orgMatrix)

    Dt = D.transpose()

    leftbottom = (0.075, 0.225)
    widthheight = (0.75, 0.75)

    fig = plt.figure(figsize=(8,8))

    vXLabels = homClusterIDs
    vYLabels = allorgs

    ax1 = fig.add_axes([0.85, leftbottom[1], 0.15, widthheight[1]])
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='right', labels=vYLabels)
    ax1.set_xticks([])
    ax1.set_yticks([])
    idx1 = Z1['leaves']

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([leftbottom[0], 0.05, widthheight[0], 0.15])
    Y = sch.linkage(Dt, method='centroid')
    Z2 = sch.dendrogram(Y, orientation='bottom', labels=vXLabels, leaf_rotation=90)
    ax2.set_xticks([])
    ax2.set_yticks([])
    idx2 = Z2['leaves']

    # Plot distance matrix.
    axmatrix = fig.add_axes([leftbottom[0], leftbottom[1], widthheight[0], widthheight[1]])

    D = D[idx1, :]
    D = D[:, idx2]

    cmap = plt.cm.get_cmap('viridis')

    cmap.set_bad(alpha=0.0)
    cmap.set_over(alpha=0.0)
    cmap.set_under(alpha=0.0)

    oRange = None
    if oRange is None:
        oRange = (np.min(D), np.max(D))

    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cmap, vmin=oRange[0], vmax=oRange[1])

    axmatrix.yaxis.set_label_position('left')
    axmatrix.xaxis.set_label_position('top')

    axmatrix.yaxis.set_ticks_position('right')
    axmatrix.xaxis.set_ticks_position('bottom')


    axmatrix.set_xticks(range(len(idx2)))
    axmatrix.set_yticks(range(len(idx1)))

    axmatrix.set_xticklabels([vXLabels[i] for i in idx2])
    axmatrix.set_yticklabels([vYLabels[i] for i in idx1])

    ax = im.axes

    for label in im.axes.xaxis.get_ticklabels():
        label.set_rotation(90)

    # Plot colorbar.
    axcolor = fig.add_axes([0.85, 0.05, 0.02, 0.15])
    minValue = np.amin(Dt)
    maxValue = np.amax(Dt)

    # cbar = mpl.colorbar.ColorbarBase(axcolor, cmap=cmap,
    #                                norm=mpl.colors.Normalize(vmin=minValue, vmax=maxValue),
    #                                orientation='horizontal')

    # pylab.colorbar(im, cax=axcolor)
    cbar = fig.colorbar(im, cax=axcolor)
    # cbar.set_label(sTitle,size=iLegendSize)

    iStep = float(cbar.vmax - cbar.vmin) / 8.0

    vTicks = []
    for i in np.arange(cbar.vmin, cbar.vmax, iStep):
        vTicks.append(i)

    if int(cbar.vmax) - vTicks[len(vTicks) - 1] < vTicks[1] - vTicks[0]:
        vTicks[len(vTicks) - 1] = cbar.vmax
    else:
        vTicks.append(cbar.vmax)

    cbar.set_ticks(vTicks)
    cbar.ax.tick_params(labelsize=8)


    plt.show()