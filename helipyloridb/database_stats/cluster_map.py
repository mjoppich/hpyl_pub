import sys,os

from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from porestat.utils.DataFrame import DataFrame

import numpy as np

from sklearn import decomposition

import matplotlib.pyplot as plt

if __name__ == '__main__':

    fileLocation = "/mnt/c/dev/data/haas/homdb/"

    homDB = HomologyDatabase.loadFromFile(fileLocation + "/combed")
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

    mc = ['4_N1-031C1', '2_N1-025A2', '14_1-20A_UB64', '13_N5-004A1', '3_N1-029C1', '11_N4-029C2', '10_N2-085C2', '1_N1-024A1']
    nmc = [x for x in allorgs if not x in mc] # and not x.startswith("15")


    allData = DataFrame()


    allData.addColumns(allorgs)
    homClusterIDs = []


    for homid in homDB.homologies:

        val = homDB.get_homology_cluster(homid)

        maxlength = 0
        for org in val:

            geneid = val[org]
            seq = genomDB.get_sequence(org, geneid)

            if len(seq) > maxlength:
                maxlength = len(seq)

        if maxlength < 80:
            continue

        allincluded = 0

        if len(val) < 2:
            continue

        for org in nmc:
            if org in val:
                allincluded += 1

        #if allincluded > 1:# len(allorgs):
        #    continue

        allincluded = 0

        for org in mc:
            if org in val:
                allincluded += 1

        #if allincluded <= 1:
        #    continue

        homClusterIDs.append(homid)

    print(allorgs)
    print(len(homClusterIDs))

    orgMatrix = []

    for org in mc + nmc:

        orgRes = []
        for homID in homClusterIDs:
            val = org in homDB.get_homology_cluster(homID)

            if val:
                orgRes.append(1)
            else:
                orgRes.append(0)

        orgMatrix.append(orgRes)

    D = np.matrix(orgMatrix)
    colors = [1 for x in mc] + [2 for x in nmc]

    # Build a forest and compute the feature importances
    forest = ExtraTreesClassifier(n_estimators=5, max_features=2, random_state=0)

    forest.fit(D, colors)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    showFeatCount = min([30, D.shape[1], len(homClusterIDs)])

    interestHomCluster = []

    for f in range(showFeatCount):
        interestHomCluster.append(homClusterIDs[f])

        print("%d. feature %d (%f) %s" % (f + 1, indices[f], importances[indices[f]], homClusterIDs[f]))


    for homID in interestHomCluster:

        val = homDB.get_homology_cluster(homID)

        mcc = sum([1 for x in val if x in mc])
        nmcc = sum([1 for x in val if x in nmc])


        print(homID, "MC", mcc, "NMC", nmcc)

    print([int(x.replace('HOMID', '')) for x in interestHomCluster])


    # Plot the feature importances of the forest
    plt.figure()
    plt.title("Feature importances")
    plt.bar(range(D.shape[1])[:showFeatCount], importances[indices][:showFeatCount],
            color="r", yerr=std[indices][:showFeatCount], align="center")
    plt.xticks(range(showFeatCount), indices[:showFeatCount])
    plt.xlim([-1, showFeatCount])
    plt.show()

    model = SelectFromModel(forest, prefit=True)
    X_new = model.transform(D)
    print(    X_new.shape)







    exit(0)


    pca = decomposition.PCA(n_components=8)
    projected = pca.fit_transform(D)

    print(projected)

    print(pca.components_)


    plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance');
    plt.show()



    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(projected[:,0], projected[:,1], projected[:, 2], c=colors, cmap=plt.cm.get_cmap('spectral', 3), alpha=0.5, marker='o')


    #for i,x in enumerate(mc+nmc):
    #    plt.annotate(x, (projected[i,0], projected[i,1], projected[i,2]))

    plt.xlabel('component 1')
    plt.ylabel('component 2')

    plt.show()

    exit(0)

    import scipy.cluster.hierarchy as sch

    # Generate random features and distance matrix.
    D = np.matrix(orgMatrix)

    Dt = D.transpose()

    leftbottom = (0.075, 0.225)
    widthheight = (0.75, 0.75)

    fig = plt.figure(figsize=(8,8))

    vXLabels = homClusterIDs
    vYLabels = mc+nmc


    idx1 = [x for x in range(0, len(vYLabels))]

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