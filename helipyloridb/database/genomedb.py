import glob
import os
from collections import defaultdict

from Bio import SeqIO


class GenomeDB:


    def __init__(self, location, loadAll = False, fileFormat='embl', fileExtension='.gb'):

        self.genomes = defaultdict(lambda: dict())
        self.location = location
        self.fileFormat = fileFormat
        self.fileExtension = fileExtension

        if loadAll:
            for file in glob.glob(location+'/*' + self.fileExtension):
                self.loadGenome(file)


    def _get_genome_id(self, file):

        gbParser = SeqIO.parse(file, self.fileFormat)
        genomeID=None

        for gb_record in gbParser:
            genomeID = gb_record.name

            if "_" in genomeID:
                genomeID = genomeID.replace('_', '')

            break

        return genomeID

    def loadGenome(self, file):


        if not os.path.isfile(file):

            if '_' in file:
                file = file.replace('_', '')

            if not os.path.isfile(file):
                loadedFile = file
                file = self.location + "/" + file + self.fileExtension

                if not os.path.isfile(file):
                    raise ValueError("Not an available genome:", file)


        if os.path.isfile(file + "e"):

            genomeID = self._get_genome_id(file)

            with open(file + "e", 'r') as infile:

                alllines = infile.readlines()

                for i in range(0, len(alllines), 2):

                    idpart = alllines[i]
                    seqpart = alllines[i+1]

                    if idpart.startswith(">"):
                        idpart = idpart[1:]

                    idpart = idpart.strip()
                    seqpart = seqpart.strip()

                    self.genomes[genomeID][idpart] = seqpart

            return



        gbParser = SeqIO.parse(file, self.fileFormat)

        with open(file + "e", 'w') as outfile:

            for gb_record in gbParser:

                genomeID = gb_record.name
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

                                                # print(nuc[frame:frame + length].translate(11))
                                                # print("%s: %s...%s - length %i, start %i, end %i, strand %i, frame %i" % (productID, longest[:30], longest[-10:], len(longest), longestStart, longestStart+len(longest)*3, strand, frame))

                                    return foundORFs

                                foundORFs = findAllOrfs(ntSeq)

                                for partProd in foundORFs:
                                    self.genomes[genomeID][partProd] = foundORFs[partProd]

                                continue

                            translation = str(ntSeq.translate())

                    if translation == None:
                        print("No Translation found for " + productID + " in genome " + genomeID)
                        continue

                    self.genomes[genomeID][productID] = translation

            for prodID in self.genomes[genomeID]:
                translation = self.genomes[genomeID][prodID]

                if translation == None:
                    print("No Translation for " + prodID + " in genome " + genomeID)
                    continue

                outfile.write(">" + str(prodID) + "\n" + translation + "\n")


            print("Loaded Genome: " + file)

    def get_sequence(self, genome, productID):
        return self.genomes.get(genome, {}).get(productID, None)

    def writeBLASTfastas(self, outpath):

        for genome in self.genomes:
            with open(outpath + "/" + genome + ".fa", 'w') as outfile:
                for seqid in self.genomes[genome]:
                    protSeq = self.genomes[genome][seqid]

                    outfile.write(">" + seqid + " " + genome + "\n")
                    outfile.write(protSeq + "\n")

    def writeCSV(self, outpath):

        with open(outpath, 'w') as outfile:

            for genome in self.genomes:
                for seqid in self.genomes[genome]:
                    protSeq = self.genomes[genome][seqid]
                    allelems = [genome, seqid, len(protSeq), protSeq]

                    allelems = [str(x).strip() for x in allelems]

                    outfile.write("\t".join(allelems) + "\n")

    def get_sequences_for_genome(self, orgname, default=None):

        if not orgname in self.genomes:
            return default

        return [x for x in self.genomes[orgname]]