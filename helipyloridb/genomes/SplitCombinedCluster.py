import argparse
import sys, os
from collections import Counter

from Bio import Phylo
from Bio.Alphabet import generic_dna
from Bio.Cluster.cluster import kmedoids
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate kmer histograms and compare for two groups', add_help=False)
    parser.add_argument('-l', '--location', type=argparse.FileType('r'), help='input', required=True)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), help='input', required=True)
    parser.add_argument('--redo', action='store_true', help='input', default=False)


    args = parser.parse_args()

    homDB = HomologyDatabase.loadFromFile(args.location.name)

    genomDB = GenomeDB(os.path.dirname(args.location.name) + "/genomes", loadAll=False)

    allorgs = homDB.get_all_organisms()

    for org in allorgs:
        genomDB.loadGenome(org)

    analyse = HomDBAnalyser(homDB, genomDB)

    maxNumberEntries = len(allorgs)

    def acceptCluster(cluster):
        if len(cluster) > maxNumberEntries:
            return True

        elif args.redo:

            clstCounter = Counter()

            for (org, orgseq) in cluster:
                clstCounter[org] += 1

            seqg2 = sum([1 for x in clstCounter if clstCounter[x] > 1])

            if seqg2 > 3:
                return True

        return False


    sizeCounter = Counter()
    for homID in homDB.homologies:
        cluster = homDB.homologies[homID]

        if acceptCluster(cluster):
            sizeCounter[len(cluster)] += 1


    cnt = 0
    for x in sizeCounter:
        cnt += sizeCounter[x]
        print(x, sizeCounter[x])

    print("Total clusters to consider", cnt)

    allOldHoms = [x for x in homDB.homologies]
    deleteHOMIDs = set()

    for homID in allOldHoms:

        cluster = homDB.homologies[homID]

        if acceptCluster(cluster):

            if args.redo:
                if not homID.startswith("sp"):
                    continue

            print(homID, len(cluster))

            #for org, seqid in cluster:
            #    genSeq = genomDB.get_sequence(org, seqid)

            #    print(">"+seqid)
            #    print(genSeq)



            aligned = analyse.cluster_align_clustalw(homID)
            for idx, cluster_ in enumerate(aligned):
                print(aligned[idx].id + "\t"+ aligned[idx].seq)


            calculator = DistanceCalculator('identity')
            constructor = DistanceTreeConstructor(calculator, 'upgma')

            tree = constructor.build_tree(aligned)
            dm = calculator.get_distance(aligned)


            tree.ladderize()  # Flip branches so deeper clades are displayed at top
            Phylo.draw_ascii(tree)

            newmat = dm.matrix
            nnmat = []
            for elem in newmat:

                nnmat.append(elem[:-1])


            newe = None
            olde = None
            nclusts = 0


            while nclusts < 8:

                nclusts += 1

                print("Testing", nclusts, "cluster(s)")
                clusterid, error, nfound = kmedoids(nnmat, nclusters=nclusts, npass=4)
                print(clusterid, error, nfound)


                olde = newe
                newe = error

                if olde != None and newe != None:

                    fact = newe/olde

                    if not args.redo:
                        if fact >= 0.8:
                            break

                    else:
                        if fact >= 0.9:
                            break

            print("Taking:")
            clusterid, error, nfound = kmedoids(nnmat, nclusters=nclusts-1, npass=4)
            print(clusterid, error, nfound)

            numClustersDetected = len(set(clusterid))


            preventDuplicate = False

            if args.redo and homID.startswith("HOMID"):
                preventDuplicate = True

            addedID = None
            if numClustersDetected > 1 and not preventDuplicate:

                deleteHOMIDs.add(homID)

                seqid2clusterentry = {}
                for elem in cluster:
                    clustid = elem[0] + "_" + elem[1]

                    seqid2clusterentry[clustid] = elem


                for clustID in set(clusterid):

                    clustSeqs = []

                    for idx, cluster in enumerate(clusterid):

                        if cluster == clustID:
                            recid = aligned[idx].id
                            recseq = str(aligned[idx].seq).replace('-', '')

                            clustSeqs.append(SeqRecord(Seq(recseq, generic_dna), id=recid, description=""))

                    clustAlign = analyse.align_clustalW_records(clustSeqs)

                    if clustAlign == None:
                        print("Empty clust align result!")
                        print(clustID)

                        for x in clustSeqs:
                            print(clustID, x)

                        continue

                    allnewClustElems = set()

                    for elem in clustAlign:
                        print(str(clustID).rjust(5), elem.id.rjust(40), elem.seq)
                        allnewClustElems.add(seqid2clusterentry[elem.id])

                    allnewClustElems = list(allnewClustElems)

                    addedID = None
                    for i in range(1, len(allnewClustElems)):

                        if addedID == None:
                            addedID = homDB.make_new_homology_relation((allnewClustElems[0], allnewClustElems[i]), None, "sp_")

                        else:
                            homDB.homologies[addedID].add(allnewClustElems[i])


                    print("Homid:", homID,"Clusters detected:", numClustersDetected, "Cluster Added", addedID)
                    print(len(homDB.homologies))
            print()
            #for idx, cluster in enumerate(clusterid):
            #    print(cluster, aligned[idx].id + "\t"+ aligned[idx].seq)


    for oldHOMID in deleteHOMIDs:
        del homDB.homologies[oldHOMID]


    homDB.save_to_file(args.output.name)