import argparse
import subprocess
import tempfile

from Bio import SeqIO, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Alphabet import generic_dna
from Bio.Cluster.cluster import kmedoids
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def align_clustalW_records(seqRecords):

    seqs = sum([1 for x in seqRecords])

    if seqs == 1:
        return seqRecords

    with tempfile.NamedTemporaryFile('w', delete=True) as tmpFastaFile, tempfile.NamedTemporaryFile('w', delete=True) as tmpMSAFile:

        try:

            # print(tmpFastaFile.name)
            # print(tmpMSAFile.name)

            SeqIO.write(seqRecords, tmpFastaFile, "fasta")
            tmpFastaFile.flush()

            clustalomega_cline = ClustalwCommandline(infile=tmpFastaFile.name, outfile=tmpMSAFile.name, output='fasta', gapopen=-1, gapext=-0.1)
            #print(clustalomega_cline)
            output = subprocess.getoutput([str(clustalomega_cline)])
            #print("msa finished", output)

            try:

                with open(tmpMSAFile.name, 'r') as fin:
                    alignment = AlignIO.read(fin, "fasta")
                    return alignment
            except:
                pass

        finally:
            pass


    return None


if __name__ == '__main__':



    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta', type=argparse.FileType('r'), required=True)

    args = parser.parse_args()


    records = [x for x in SeqIO.parse(args.fasta, 'fasta')]


    aligned = align_clustalW_records(records)

    seqCount = 0

    for idx, cluster in enumerate(aligned):
        seqCount += 1
        #print(aligned[idx].id + "\t" + aligned[idx].seq)

    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'upgma')

    tree = constructor.build_tree(aligned)
    dm = calculator.get_distance(aligned)

    #tree.ladderize()  # Flip branches so deeper clades are displayed at top
    #Phylo.draw_ascii(tree)

    newmat = dm.matrix
    nnmat = []
    for elem in newmat:
        nnmat.append(elem[:-1])

    newe = None
    olde = None

    nclusts = 0

    clustSize2IDs = {}

    while nclusts < 8 and nclusts < seqCount:

        nclusts += 1

        #print("Testing", nclusts, "cluster(s)")
        clusterid, error, nfound = kmedoids(nnmat, nclusters=nclusts, npass=4)
        #print(clusterid, error, nfound)

        clustSize2IDs[nclusts] = (clusterid, error)

        olde = newe
        newe = error

        if olde != None and newe != None:

            fact = newe / olde

            if fact >= 0.8:
                nclusts -= 1
                break

    #print("Taking", nclusts - 1, "clusters")
    clusterid, error = clustSize2IDs[nclusts]
    #print(clusterid, error)

    for clustID in set(clusterid):

        clustSeqs = []

        for idx, cluster in enumerate(clusterid):

            if cluster == clustID:
                recid = aligned[idx].id
                recseq = str(aligned[idx].seq).replace('-', '')

                clustSeqs.append(SeqRecord(Seq(recseq, generic_dna), id=recid, description=""))

        clustAlign = align_clustalW_records(clustSeqs)

        for elem in clustAlign:
            print(str(clustID).rjust(5), elem.id.rjust(30), elem.seq)