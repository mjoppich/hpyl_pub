import subprocess

import os
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline, TCoffeeCommandline, ClustalwCommandline
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase

import tempfile


class HomDBAnalyser:

    def __init__(self, homDB, genomDB):

        assert(isinstance(homDB, HomologyDatabase))
        self.homDB = homDB

        assert(isinstance(genomDB, GenomeDB))
        self.genomDB = genomDB

        for org in self.homDB.get_all_organisms():
            self.genomDB.loadGenome(org)

    def cluster_align(self, clusterID, allowedOrganisms=None):

        alignSeqs = []
        homCluster = self.homDB.get_cluster(clusterID)

        for org in homCluster:
            if allowedOrganisms == None or org in allowedOrganisms:
                for seqid in homCluster[org]:
                    alignSeqs.append((org, seqid))

        seqRecords = []
        seqID2Element = {}
        for org, seqid in alignSeqs:
            genSeq = self.genomDB.get_sequence(org, seqid)

            seqRecID = "_".join([org, seqid])

            seqID2Element[seqRecID] = self.genomDB.get_element(org, seqid)
            seq = SeqRecord(Seq(genSeq, generic_dna), id=seqRecID, description="")

            seqRecords.append(seq)


        with tempfile.NamedTemporaryFile('w', delete=True) as tmpFastaFile, tempfile.NamedTemporaryFile('w', delete=True) as tmpMSAFile:

            try:

                #print(tmpFastaFile.name)
                #print(tmpMSAFile.name)

                SeqIO.write(seqRecords, tmpFastaFile, "fasta")
                tmpFastaFile.flush()

                clustalomega_cline = ClustalOmegaCommandline(infile=tmpFastaFile.name, outfile=tmpMSAFile.name, force=True, outfmt='fa',
                                                             verbose=True, auto=True)

                clustalomega_cline = str(clustalomega_cline)
                clustalomega_cline += " --full --distmat-out /tmp/clustalodm"
                print(clustalomega_cline)
                output = subprocess.getoutput([str(clustalomega_cline)])
                print("Clustalomega finished")

                with open(tmpMSAFile.name, 'r') as fin:
                    alignment = AlignIO.read(fin, "fasta")
                    return alignment

            finally:
                pass


        return None


    def cluster_align_clustalw(self, clusterID, allowedOrganisms=None):

        alignSeqs = []
        homCluster = self.homDB.get_cluster(clusterID)

        for org in homCluster:
            if allowedOrganisms == None or org in allowedOrganisms:
                for seqid in homCluster[org]:
                    alignSeqs.append((org, seqid))

        seqRecords = []
        seqID2Element = {}
        for org, seqid in alignSeqs:
            genSeq = self.genomDB.get_sequence(org, seqid)

            seqRecID = "_".join([org, seqid])

            seqID2Element[seqRecID] = self.genomDB.get_element(org, seqid)
            seq = SeqRecord(Seq(genSeq, generic_dna), id=seqRecID, description="")

            seqRecords.append(seq)


        with tempfile.NamedTemporaryFile('w', delete=True) as tmpFastaFile, tempfile.NamedTemporaryFile('w', delete=True) as tmpMSAFile:

            try:

                #print(tmpFastaFile.name)
                #print(tmpMSAFile.name)

                SeqIO.write(seqRecords, tmpFastaFile, "fasta")
                tmpFastaFile.flush()

                clustalomega_cline = ClustalwCommandline(infile=tmpFastaFile.name, outfile=tmpMSAFile.name, output='fasta', gapopen=-10, gapext=-0.01)
                print(clustalomega_cline)
                output = subprocess.getoutput([str(clustalomega_cline)])
                print("msa finished")

                try:

                    with open(tmpMSAFile.name, 'r') as fin:
                        alignment = AlignIO.read(fin, "fasta")
                        return alignment
                except:
                    pass

            finally:
                pass


        return None


