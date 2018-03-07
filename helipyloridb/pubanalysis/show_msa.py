import glob
from Bio import AlignIO
from vis.msa_viewer import MSAViewer

for file in glob.glob("/mnt/c/ownCloud/data/hpyloriDB/luisa_msa/*.aln"):

    print(file)

    align = AlignIO.read(file, "clustal")

    allRecords = []

    for record in align:
        allRecords.append( (record.id, record.seq) )


    partLen = 1024

    recLen = len(allRecords[0][1])
    if recLen > partLen:

        replRecs = {}

        for i in range(0, recLen, partLen):
            allpartSeqs = []

            for (recid, seq) in allRecords:
                allpartSeqs.append((recid, seq[i:i+partLen]))

            replRecs["Part" + str(i)] = allpartSeqs

        allRecords = replRecs



    MSAViewer.makeMSA(allRecords, file + ".html")
