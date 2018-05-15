import glob

import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase

if __name__ == '__main__':
    baseDIR = '/mnt/c/dev/data/haas/homdb/'

    homDB = HomologyDatabase.loadFromFile(baseDIR + "/combed")

    homs = ['HOMID1448', 'HOMID1742', 'HOMID1692', 'HOMID1795', 'HOMID2024', 'HOMID2027', 'HOMID1621', 'HOMID1338', 'HOMID1672', 'HOMID1693']
    homs = ['HOMID1286']
    homs = ['HOMID933', 'HOMID1354', 'HOMID1792', 'HOMID1621', 'HOMID1165', 'HOMID2171', 'HOMID283']
    homs = ['HOMID933', 'HOMID1354']

    promLen = 35

    for homid in homs:

        alignSeqs = []
        homCluster = homDB.get_cluster(homid)

        for x in homCluster:
            print(x, homCluster[x])


        print()
        print()

        for file in glob.glob(baseDIR + "/genomes/*.gb"):

            org = os.path.basename(file).replace('.gb', '')

            orgLT = homCluster.get(org, None)

            if orgLT == None:
                continue

            for rec in SeqIO.parse(file, 'embl'):

                for feature in rec.features:

                    if not feature.type == 'CDS':
                        continue

                    fltContained = any([x in orgLT for x in feature.qualifiers['locus_tag']])

                    if fltContained:

                        #print(feature)
                        if feature.location.strand == 1:
                            print(rec.seq[feature.location.start-promLen:feature.location.start], org, feature.qualifiers['locus_tag'], homid)
                        else:
                            rawSeq = rec.seq[feature.location.end:feature.location.end+promLen]
                            rawSeq = rawSeq.reverse_complement()
                            print(rawSeq, org, feature.qualifiers['locus_tag'], homid)
