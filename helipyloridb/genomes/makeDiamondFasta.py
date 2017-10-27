import glob

from Bio import SeqIO

from utils import fileLocation

for file in glob.glob(fileLocation + '/genomes/*.gb'):

    gbParser = SeqIO.parse(file, "embl")

    for gb_record in gbParser:

        genomeID = gb_record.name
        allSeqs = {}

        mainFeature = None

        for feature in gb_record.features:

            if feature.type.upper() == 'SOURCE':
                mainFeature = feature

            if not feature.type.upper() in ['CDS', 'GENE']:
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

                if feature.type.upper() == 'GENE' and mainFeature != None:
                    ntSeq = feature.location.extract(gb_record).seq

                    modLen = len(ntSeq) % 3

                    if modLen != 0:

                        def findAllOrfs(seq):

                            foundORFs = {}

                            for strand, nuc in [(+1, seq)]:
                                for frame in range(3):
                                    length = 3 * ((len(seq) - frame) // 3)  # Multiple of three
                                    longest = ""
                                    longestStart = 0

                                    startpos = frame

                                    for pro in nuc[frame:frame + length].translate(11).split("*"):

                                        if len(pro) > len(longest):
                                            longest = pro
                                            longestStart = startpos

                                        startpos += len(pro) * 3 + 1

                                    if (len(longest) > 15) and (frame == 0 or frame == modLen):

                                        foundORFs[productID + "_" + str(frame)] = str(longest)

                                        #print(nuc[frame:frame + length].translate(11))
                                        #print("%s: %s...%s - length %i, start %i, end %i, strand %i, frame %i" % (productID, longest[:30], longest[-10:], len(longest), longestStart, longestStart+len(longest)*3, strand, frame))

                            return foundORFs

                        foundORFs = findAllOrfs(ntSeq)

                        for partProd in foundORFs:
                            allSeqs[partProd] = foundORFs[partProd]

                        continue

                    translation = str(ntSeq.translate())
                else:
                    continue


            allSeqs[productID] = translation

        with open(fileLocation+'/genomes/' + genomeID + '.fa', 'w') as outfile:

            for prod in sorted([x for x in allSeqs]):
                seq = allSeqs[prod]
                outfile.write(">" + prod + " " + genomeID + "\n")
                outfile.write(seq + "\n")

            print("Done: " + genomeID)
