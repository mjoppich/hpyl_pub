import editdistance

from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase

if __name__ == '__main__':

    genomeLocation = '/home/users/joppich/ownCloud/data/hpyloriDB/genomes/'

    homDB = HomologyDatabase.loadFromFile("/home/proj/projekte/dataintegration/hpyloriDB/hpp12.homdb")
    genDB = GenomeDB(genomeLocation)

    for homGroup in homDB.homologies:

        entries = homDB.homologies[homGroup]

        allSeqs = []

        for seqID in entries:

            if not seqID[0] in genDB.genomes:
                genDB.loadGenome(genomeLocation + "/" + seqID[0] + ".gb")

            seq = genDB.get_sequence(seqID[0], seqID[1])

            allSeqs.append(seq)

        if len(allSeqs) == 0:
            continue

        seqEnds = []

        for seq in allSeqs:

            startSeq = max(len(seq)-10, 0)
            endSeq = len(seq)

            seqEnds.append( seq[startSeq:endSeq] )

        if len(set(seqEnds)) > 1:

            listSeqs = list(seqEnds)
            distances = []

            for i in range(0, len(seqEnds)):
                for j in range(i+1, len(seqEnds)):

                    seqI = seqEnds[i]
                    seqJ = seqEnds[j]

                    dist = editdistance.eval(seqI, seqJ)
                    distances.append( (seqI, seqJ, dist) )

            print(homGroup + "\t" + str(len(allSeqs)) + "\t" + str(len(set(seqEnds))) + "\t" + ",".join(set(seqEnds)) + "\t" + ",".join([str(x) for x in distances]))
