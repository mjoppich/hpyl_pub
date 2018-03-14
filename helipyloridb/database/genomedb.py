import glob
import json
import os
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

from utils import fileLocation

class GenomeDBEntry:

    def __init__(self, organismID, organismName, recordID, entryID, seqAA, seqNT, genomicStart, genomicEnd, genomicStrand, entryNames):
        self.organismID = organismID
        self.organismName = organismName
        self.recordID = recordID
        self.entryID = entryID
        self.seqAA = str(seqAA)
        self.seqNT = str(seqNT)
        self.genomicStart = int(genomicStart)
        self.genomicEnd = int(genomicEnd)
        self.genomicStrand = genomicStrand
        self.entryNames = entryNames

    def __str__(self):

        return "\t".join( [str(x) for x in [ self.organismID,
                            self.organismName,
                            self.recordID,
                            self.entryID,
                            self.seqAA,
                            self.seqNT,
                            self.genomicStart,
                            self.genomicEnd,
                            self.genomicStrand,
                            self.entryNames
                            ]] )

    def toJSON(self):

        retDict = {}

        for x in self.__dict__:
            retDict[x] = self.__dict__[x]

        return retDict

    @classmethod
    def fromLine(cls, line):

        aline = line.strip().split('\t')
        idx2content = defaultdict(lambda: None)

        for idx, con in enumerate(aline):
            idx2content[idx] = con

        return GenomeDBEntry(
            idx2content[0],
            idx2content[1],
            idx2content[2],
            idx2content[3],
            idx2content[4],
            idx2content[5],
            idx2content[6],
            idx2content[7],
            idx2content[8],
            eval(idx2content[9])
        )


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

    def loadGenome(self, file, force=False):


        if not os.path.isfile(file):

            if '_' in file:
                file = file.replace('_', '')

            if not os.path.isfile(file):
                loadedFile = file
                file = self.location + "/" + file + self.fileExtension

                if not os.path.isfile(file):
                    raise ValueError("Not an available genome:", file)


        if force==False and os.path.isfile(file + "e"):

            genomeID = self._get_genome_id(file)

            with open(file + "e", 'r') as infile:

                for line in infile:
                    genomeEntry = GenomeDBEntry.fromLine(line)
                    self.genomes[genomeID][genomeEntry.entryID] = genomeEntry

            return

        def makeEntryNames(feature):

            allnotes = feature.qualifiers.get('note', [])

            noteInfo = set()
            for x in allnotes:
                notes = x.split('; ')

                for note in notes:
                    noteInfo.add(note)

            noteInfo = list(noteInfo)

            return [x for x in [feature.qualifiers.get('locus_tag', [None])[0], feature.qualifiers.get('product', [None])[0]] + noteInfo if x != None]

        gbParser = SeqIO.parse(file, self.fileFormat)

        with open(file + "e", 'w') as outfile:

            for gb_record in gbParser:

                genomeID = gb_record.name
                mainFeature = None

                for feature in gb_record.features:

                    if feature.type.upper() == 'SOURCE':
                        mainFeature = feature

                    if not feature.type.upper() in ['GENE', 'CDS']:
                        continue

                    locTag = feature.qualifiers.get('locus_tag', [None])[0]
                    proID = feature.qualifiers.get('protein_id', [None])[0]
                    translation = feature.qualifiers.get('translation', [None])[0]
                    ntSeq = feature.location.extract(gb_record).seq

                    productID = locTag if locTag != None else proID

                    if productID == None:
                        print("CDS without id:")
                        print(feature)
                        continue

                    if translation == None:

                        if feature.type.upper() == 'GENE' and mainFeature != None:
                            ntSeq = feature.location.extract(gb_record).seq # TODO is this already RC on - ?

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

                                                startpos += len(pro) * 3 + 3

                                            if (len(longest) > 15) and (frame == 0 or frame == modLen):
                                                testSeq = nuc[longestStart:longestStart+3*len(str(longest))]
                                                foundORFs[productID + "_" + str(frame)] = (str(longest), longestStart, longestStart + 3*len(str(longest)))

                                                # print(nuc[frame:frame + length].translate(11))
                                                # print("%s: %s...%s - length %i, start %i, end %i, strand %i, frame %i" % (productID, longest[:30], longest[-10:], len(longest), longestStart, longestStart+len(longest)*3, strand, frame))

                                    return foundORFs

                                foundORFs = findAllOrfs(ntSeq)

                                for partProd in foundORFs:

                                    partProdData = foundORFs[partProd]

                                    gstart = feature.location.start

                                    if feature.location.strand == '-1':
                                        gstart -= partProdData[1]
                                        gend = gstart- partProdData[2]

                                    else:
                                        gstart = feature.location.start + partProdData[1]
                                        gend = gstart + partProdData[2]


                                    subNtSeq = ntSeq[gstart:gend]

                                    partEntry = GenomeDBEntry(
                                        organismID=genomeID,
                                        organismName=", ".join(mainFeature.qualifiers['organism']),
                                        recordID=gb_record.id,
                                        entryID=partProd,
                                        seqAA=partProdData[0],
                                        seqNT=subNtSeq,
                                        genomicStart=gstart,
                                        genomicEnd=gend,
                                        genomicStrand=feature.location.strand,
                                        entryNames=makeEntryNames(feature),
                                    )

                                    self.genomes[genomeID][partProd] = partEntry

                                continue

                            translation = str(ntSeq.translate())
                    else:
                        translation += "*"

                    if translation == None:
                        print("No Translation found for " + productID + " in genome " + genomeID)
                        continue

                    transAA = Seq(str(ntSeq)).translate()

                    if str(transAA) != translation:
                        print("error in translation")
                        print(productID)
                        print(str(transAA))
                        print(translation)
                        print()

                    entry = GenomeDBEntry(
                        organismID=genomeID,
                        organismName=", ".join(mainFeature.qualifiers['organism']),
                        recordID=gb_record.id,
                        entryID=productID,
                        entryNames=makeEntryNames(feature),
                        seqAA=translation,
                        seqNT=ntSeq,
                        genomicStart=feature.location.start,
                        genomicEnd=feature.location.end,
                        genomicStrand=feature.location.strand
                    )


                    self.genomes[genomeID][productID] = entry

            for prodID in self.genomes[genomeID]:
                genomicEntry = self.genomes[genomeID][prodID]

                if genomicEntry == None:
                    print("No Translation for " + prodID + " in genome " + genomeID)
                    continue

                outfile.write( str(genomicEntry) +"\n" )


            print("Loaded Genome: " + file)

    def get_sequence(self, genome, productID):
        genElem = self.genomes.get(genome, {}).get(productID, None)

        return None if genElem == None else genElem.seqAA

    def writeBLASTfastas(self, outpath):

        for genome in self.genomes:
            with open(outpath + "/" + genome + ".fa", 'w') as outfile:
                for seqid in self.genomes[genome]:
                    protSeq = self.genomes[genome][seqid].seqAA

                    outfile.write(">" + seqid + " " + genome + "\n")
                    outfile.write(protSeq + "\n")

    def writeCSV(self, outpath):

        with open(outpath, 'w') as outfile:

            for genome in self.genomes:
                for seqid in self.genomes[genome]:
                    protSeq = self.genomes[genome][seqid].seqAA
                    allelems = [genome, seqid, len(protSeq), protSeq]

                    allelems = [str(x).strip() for x in allelems]

                    outfile.write("\t".join(allelems) + "\n")

    def get_sequences_for_genome(self, orgname, default=None):

        if not orgname in self.genomes:
            return default

        return [x for x in self.genomes[orgname]]

    def get_element(self, genome, productID):
        genElem = self.genomes.get(genome, {}).get(productID, None)

        return None if genElem == None else genElem

if __name__ == '__main__':

    genomDB = GenomeDB(fileLocation + "/genomes/")
    genomDB.loadGenome('AE000511', force=True)

    res = genomDB.get_element('AE000511', 'HP_1474')

    print(res)
    print(res.toJSON())

    print(json.dumps(res.toJSON(), indent=4, sort_keys=True))