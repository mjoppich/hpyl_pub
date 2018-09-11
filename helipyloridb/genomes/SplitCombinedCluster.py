import argparse
import random
import sys, os
from collections import Counter
from datetime import datetime

from Bio import Phylo
from Bio.Alphabet import generic_dna
from Bio.Cluster.cluster import kmedoids
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase
import statistics

from utils.Parallel import MapReduce

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from analysis.homologybuilder import HomologyBuilder


import numpy

def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.

    A cell (i,j) in the array is the length of the branch between allclades[i]
    and allclades[j], if a branch exists, otherwise infinity.

    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a NumPy 2D array.
    """
    allclades = list(tree.find_clades(order='level'))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    distmat = numpy.repeat(numpy.inf, len(allclades)**2)
    distmat.shape = (len(allclades), len(allclades))
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            if child.branch_length:
                distmat[lookup[parent], lookup[child]] = child.branch_length
    if not tree.rooted:
        distmat += distmat.transpose
    return (allclades, numpy.matrix(distmat))

def makeClusterID(clusterID, existingIDs):

    new_clust_name = "cluster_" + str(clusterID)
    while new_clust_name in existingIDs:
        new_clust_name = "cluster_" + str(random.randrange(1, 1000))

    return new_clust_name


def processClusterInitial(homID):

    seqid2orgtuple, seqs = analyse.get_cluster_records(homID, allowedOrganisms=restrictOrgs)

    return seqid2orgtuple, processCluster(seqs)



def processCluster(seqs, refine=True, startClusterCount=1):


    aligned = analyse.align_clustalW_records(seqs)

    if aligned == None:
        return None

    for idx, cluster_ in enumerate(aligned):
        print(aligned[idx].id.rjust(30) + "\t" + aligned[idx].seq)


    if len(seqs) == 1:
        print("Single Sequence")
        return None


    calculator = DistanceCalculator('blosum80')
    dm = calculator.get_distance(aligned)


    newmat = dm.matrix
    nnmat = []
    for elem in newmat:
        nnmat.append(elem[:-1])

    newe = None
    olde = None
    nclusts = startClusterCount-1

    clustSize2IDs = []

    clustPasses = 6

    if len(seqs) < 10:
        clustPasses = 12

    breakForLowDissimBetweenClusters = False

    while nclusts < len(seqs):

        nclusts += 1

        print("Testing", nclusts, "cluster(s)")
        clusterid, error, nfound = kmedoids(nnmat, nclusters=nclusts, npass=clustPasses)
        print(clusterid, error, nfound)

        clustSize2IDs.append((nclusts, clusterid, error))



        if error == 0.0:
            # overfit?
            if refine==True:
                nclusts += 1

            break


        olde = newe
        newe = error

        canBreak = False

        if olde != None and newe != None:

            fact = newe / olde

            if olde < 10:

                if newe < maxAllowedDissimWithinCluster:
                    canBreak = True
                    breakForLowDissimBetweenClusters = True
                    break
                else:
                    if not args.redo:
                        if fact >= 0.8:
                            canBreak = True


                    else:
                        if fact >= 0.9:
                            canBreak = True

            else:

                if abs(olde - newe) <= 1.0:
                    canBreak = True


        maxDiffInCluster = 0.0

        for clustID in set(clusterid):
            clustSeqs = []
            for idx, subcluster in enumerate(clusterid):

                if subcluster == clustID:
                    clustSeqs.append(aligned[idx].id)


            for i in range(0, len(clustSeqs)):
                for j in range(i, len(clustSeqs)):
                    idi = clustSeqs[i]
                    idj = clustSeqs[j]

                    diffScore = dm[(idi, idj)]

                    maxDiffInCluster = max([maxDiffInCluster, diffScore])

        if maxDiffInCluster < maxAllowedDissimWithinCluster:
            canBreak = True

        if canBreak:

            if refine==True and nclusts < 3:
                continue
            break

    allClusters = {}
    reprocessClusters = []

    if refine:
        print("In Refine Mode")

    if len(clustSize2IDs) > 1:
        nclusts, clusterid, error = clustSize2IDs[-2]

        if nclusts == 1 and refine:
            print("Saving from endless loop")
            nclusts, clusterid, error = clustSize2IDs[-1]

    else:
        nclusts, clusterid, error = clustSize2IDs[-1]

    print("Taking", nclusts, "clusters")
    print(clusterid, error)

    if nclusts == len(seqs):
        print("all seqs have own cluster ...", )

    for clustID in set(clusterid):

        clustSeqs = []

        for idx, subcluster in enumerate(clusterid):

            if subcluster == clustID:
                recid = aligned[idx].id
                recseq = str(aligned[idx].seq).replace('-', '')

                clustSeqs.append(SeqRecord(Seq(recseq, generic_dna), id=recid, description=""))

        """
        find a unique cluster id !
        """
        new_clust_name = makeClusterID(clustID, [x for x in allClusters])

        allClusters[new_clust_name] = clustSeqs

        if len(clustSeqs) > 1:

            if nclusts > 1:
                print("calling clustalo for", len(clustSeqs), "seqs")
                clustAlign = analyse.align_clustalo_records(clustSeqs)

                for elem in clustAlign:
                    print(str(new_clust_name).rjust(10), str(clustID).rjust(5), elem.id.rjust(40), elem.seq)

                dm = calculator.get_distance(clustAlign)


            allPWDists = sorted([y for x in dm.matrix for y in x])
            maxDissimWithinCluster = max(allPWDists)

            if breakForLowDissimBetweenClusters:
                maxDissimWithinCluster = statistics.median(allPWDists)

            print("Max Dissim Detected ", "median" if breakForLowDissimBetweenClusters else "", maxDissimWithinCluster)
            for ridx, row in enumerate(dm.matrix):
                for lidx, elem in enumerate(row):
                    if elem >= maxDissimWithinCluster:
                        print(dm.names[ridx], dm.names[lidx], elem)

            print()
            print()

            maxAllowedDissim = maxAllowedDissimByLength.get(len(clustSeqs), maxAllowedDissimWithinCluster)

            if maxDissimWithinCluster > maxAllowedDissim:


                reprocessClusters.append(new_clust_name)
                print("Max dissim within cluster", maxDissimWithinCluster, ">", maxAllowedDissim, new_clust_name, [seq.id for seq in clustSeqs])


    for clustID in reprocessClusters:

        clustSeqs = allClusters[clustID]

        startClusters = 2
        if len(clustSeqs) < 4:
            startClusters = 1

        newClusters = processCluster(allClusters[clustID], refine=True, startClusterCount=startClusters)

        if newClusters != None:

            del allClusters[clustID]

            for cID in newClusters:

                newcid = makeClusterID(clustID, [x for x in allClusters])
                allClusters[newcid] = newClusters[cID]

    return allClusters

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=argparse.FileType('r'), help='input', required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='input', required=True)
    parser.add_argument('--redo', action='store_true', help='input', default=False)


    restrictOrgs = ['AE001439', 'AE000511', 'CP001217']
    restrictOrgs = None


    args = parser.parse_args()

    print("Loading Hom DB")
    homDB = HomologyDatabase.loadFromFile(args.location.name)

    print("Loading Genomes")
    genomDB = GenomeDB(os.path.dirname(args.location.name) + "/genomes", loadAll=False)

    allorgs = homDB.get_all_organisms()

    if restrictOrgs:
        allorgs = restrictOrgs

    for org in allorgs:
        genomDB.loadGenome(org)


    print("Loading HomDB analyser")
    analyse = HomDBAnalyser(homDB, genomDB, loadAll=False)

    maxNumberEntries = len(allorgs)

    maxAllowedDissimWithinCluster = 0.25

    maxAllowedDissimByLength = {}
    maxAllowedDissimByLength[5] = 0.3
    maxAllowedDissimByLength[4] = 0.4
    maxAllowedDissimByLength[3] = 0.45
    maxAllowedDissimByLength[2] = 0.5


    changeDB = False


    def acceptCluster(cluster):

        allEntryCount = [len(cluster[x]) for x in cluster if x in allorgs]

        if len(allEntryCount) == 0:
            return False

        if max(allEntryCount) > 2:
            return True

        if sum(allEntryCount) > len(allorgs):
            return True

        elif args.redo:

            clstCounter = Counter()

            for (org, orgseq) in cluster:
                clstCounter[org] += 1

            seqg2 = sum([1 for x in clstCounter if clstCounter[x] > 1])

            if seqg2 > 3:
                return True

        return False


    print("Scanning clusters")

    acceptedClusters = []
    sizeCounter = Counter()
    for homID in homDB.homologies:
        cluster = homDB.get_cluster(homID)

        if acceptCluster(cluster):

            acceptedClusters.append(homID)
            sizeCounter[len(cluster)] += 1


    cnt = 0
    for x in sizeCounter:
        cnt += sizeCounter[x]
        print(x, sizeCounter[x])



    print("Total clusters to consider", cnt)
    print("Clusters for reprocessing", len(acceptedClusters))
    print(acceptedClusters)

    with open(args.location.name + ".progress", "w") as myfile:
        myfile.write(str(datetime.now()) + " Starting SplitCombinedCluster\n")


    def logToFile(message):

        with open(args.location.name + ".progress", "a") as myfile:
            myfile.write(str(datetime.now()) + " " + message + "\n")


    def procFunc(datas, env):

        res = []

        for homID in datas:

            logToFile("Processing Cluster " + homID)

            cluster = homDB.get_cluster(homID)

            if acceptCluster(cluster):

                if restrictOrgs:
                    clusterold = cluster
                    cluster = {}
                    for x in restrictOrgs:
                        if x in clusterold:
                            cluster[x] = clusterold[x]

                if args.redo:
                    if not homID.startswith("sp"):
                        continue

                print(homID, len(cluster))

                # for org, seqid in cluster:
                #    genSeq = genomDB.get_sequence(org, seqid)

                #    print(">"+seqid)
                #    print(genSeq)

                """
    
                returns a map: clusterID => sequence
    
                """
                seqid2orgtuple, allClusterSeqs = processClusterInitial(homID)

                res.append(  (homID, seqid2orgtuple, allClusterSeqs) )

        return res


    def joinFunc(oldres, newres, env):

        for res in newres:

            (homID, seqid2orgtuple, allClusterSeqs) = res

            logToFile("Finished Clustering for hom ID " + str(homID) + ": " + str(len(allClusterSeqs)))

            seqCount = 0
            clusterSizes = Counter()

            for clusterID in allClusterSeqs:

                print("Finish clustalo for", len(allClusterSeqs[clusterID]), "seqs")
                clustAlign = analyse.align_clustalo_records(allClusterSeqs[clusterID])

                clusterSizes[len(clustAlign)] += 1

                for elem in clustAlign:
                    print(str(clusterID).rjust(20), elem.id.rjust(40), elem.seq)

                    seqCount += 1

            logToFile("Finished Clustering for hom ID " + str(homID) + " with " + str(seqCount) + " sequences : " + str(
                clusterSizes))

            preventDuplicate = False
            if args.redo and homID.startswith("HOMID"):
                preventDuplicate = True

            addedID = None
            if len(allClusterSeqs) > 1 and not preventDuplicate:

                for clusterID in allClusterSeqs:

                    allnewClustElems = set()
                    for seq in allClusterSeqs[clusterID]:
                        allnewClustElems.add(seqid2orgtuple[seq.id])

                    allnewClustElems = list(allnewClustElems)

                    if changeDB:

                        addedID = None
                        for i in range(1, len(allnewClustElems)):

                            if addedID == None:
                                addedID = homDB.make_new_homology_relation((allnewClustElems[0], allnewClustElems[i]),
                                                                           {'original_cluster': homID}, "sp_")

                                logToFile("##" + str(homID) + "\t" + str(addedID) + "##")


                            else:
                                homDB.homologies[addedID].add(allnewClustElems[i])

                        print("Homid:", homID, "Clusters detected:", len(allClusterSeqs), "Cluster Added", addedID)
                        print(len(homDB.homologies))
                        print()

                if changeDB:
                    del homDB.homologies[homID]
            # for idx, cluster in enumerate(clusterid):
            #    print(cluster, aligned[idx].id + "\t"+ aligned[idx].seq)


        homDB.save_to_file(args.output.name)


    logToFile("Processing Clusters ("+str(len(acceptedClusters))+") " + str(acceptedClusters))

    ll = MapReduce(4)
    ll.exec(acceptedClusters, procFunc, None, 1, joinFunc)