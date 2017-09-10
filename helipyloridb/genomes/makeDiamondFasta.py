import glob

from Bio import SeqIO

for file in glob.glob('../../../genomes/*.gb'):

    gbParser = SeqIO.parse(file, "embl")

    for gb_record in gbParser:

        genomeID = gb_record.name
        allSeqs = {}

        for feature in gb_record.features:

            if not feature.type.upper() == 'CDS':
                continue

            locTag = feature.qualifiers.get('locus_tag', [None])[0]
            proID = feature.qualifiers.get('protein_id', [None])[0]
            translation = feature.qualifiers.get('translation', [None])[0]

            productID = locTag if locTag != None else proID

            if productID == None:
                print("CDS without id:")
                print(feature)
                continue

            if translation == None:
                continue

            allSeqs[productID] = translation

        with open('../../../genomes/' + genomeID + '.fa', 'w') as outfile:

            for prod in allSeqs:
                seq = allSeqs[prod]
                outfile.write(">" + prod + " " + genomeID + "\n")
                outfile.write(seq + "\n")

            print("Done: " + genomeID)
