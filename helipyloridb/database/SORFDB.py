import re
from collections import defaultdict

from openpyxl import load_workbook


class SORFDB:

    def __init__(self):

        self.sorfs = defaultdict(list)


    def find_gene(self, geneID, defaultval=list()):

        retIDs = []

        for x in self.sorfs:

            if geneID == self.sorfs[x]['SORF_ASSOC_LT']:
                retIDs.append(x)

        return defaultval if len(retIDs) == 0 else retIDs

    def get_sorf_infos(self, sorfID):

        if not sorfID in self.sorfs:
            return None

        sorfBase = self.sorfs[sorfID]

        return sorfBase


    @classmethod
    def from_cs_sorfs(cls, filelocation="/mnt/c/ownCloud/data/hpyloriDB/sharma/sorfs.xlsx"):

        sorfDB = SORFDB()

        wb = load_workbook(filelocation, data_only=True)
        ws = wb.active

        curSORFID = 0

        sorfLTNameCol = 0
        sorfStartCol = 1
        sorfStrandCol = 2

        sorfAltNameCol = 3
        sorfAntisenseLTCol = 4
        sorfAssociatedLTCol = 5

        sorfNameCol = 6
        sorfEndCol = 7

        sorfTransTermCol = 8
        sorfTransTermEndCol = 9

        sorfMax3EndCol = 10
        sorfSRNAClassCol = 11
        sorfPrimaryTSSCol = 12
        sorfSecondaryTSSCol = 13
        sorfInternalTSSCol = 14

        sorfAntisenseCol = 15
        sorfPutativeSRNACol = 16
        sorfPutativeASRNACol = 17
        sorfMaybe5UTRCol = 18
        sorfSeqCol = 19

        def transform_lt(inLT):

            outLT = inLT
            if inLT != None and inLT[0:2] == 'HP':
                outLT = inLT[0:2] + "_" + inLT[2:]

            return outLT


        for row in ws:
            if row[0].row < 9:
                continue


            sorfLT = row[sorfLTNameCol].value
            sorfStart = row[sorfStartCol].value
            sorfEnd = row[sorfEndCol].value

            sorfStrand = row[sorfStrandCol].value
            sorfAltName = row[sorfAltNameCol].value
            sorfAntisenseLT = row[sorfAntisenseLTCol].value
            sorfAssociatedLT = row[sorfAssociatedLTCol].value
            sorfName = row[sorfNameCol].value

            sorfTransTerm = row[sorfTransTermCol].value
            sorfTransTermEnd = row[sorfTransTermEndCol].value

            sorfMax3End = row[sorfMax3EndCol].value
            sorfSeq = row[sorfSeqCol].value

            sorfLT = transform_lt(sorfLT)
            sorfAntisenseLT = transform_lt(sorfAntisenseLT)
            sorfAssociatedLT = transform_lt(sorfAssociatedLT)

            sorfInfos = {
                'LOCUS_TAG': sorfLT,
                'SORF_START': sorfStart,
                'SORF_END': sorfEnd,
                'SORF_STRAND': sorfStrand,
                'SORF_ALTNAME': sorfAltName,
                'SORF_ANTISENSE_LT': sorfAntisenseLT,
                'SORF_ASSOC_LT': sorfAssociatedLT,
                'SORF_NAME': sorfName,
                'SORF_TRANSTERM_TERM': sorfTransTerm,
                'SORF_TRANSTERM_END': sorfTransTermEnd,
                'SORF_MAX_3END': sorfMax3End,
                'SORF_SEQ': sorfSeq
            }


            sorfProps = []


            if row[sorfSRNAClassCol].value == 1:
                sorfProps.append('class')

            if row[sorfPrimaryTSSCol].value == 1:
                sorfProps.append('primary')

            if row[sorfSecondaryTSSCol].value == 1:
                sorfProps.append('secondary')

            if row[sorfInternalTSSCol].value == 1:
                sorfProps.append('internal')

            if row[sorfAntisenseCol].value == 1:
                sorfProps.append('antisense')

            if row[sorfPutativeSRNACol].value == 1:
                sorfProps.append('putative sRNA')

            if row[sorfPutativeASRNACol].value == 1:
                sorfProps.append('putative asRNA')

            if row[sorfMaybe5UTRCol].value == 1:
                sorfProps.append('maybe 5\' UTR')


            curSORF = 'SORF' + str(curSORFID)
            curSORFID += 1

            sorfInfos['PROPS'] = sorfProps


            sorfDB.sorfs[curSORF] = sorfInfos

        return sorfDB


if __name__ == '__main__':

    sorfDB = SORFDB.from_cs_sorfs()

    sorfIDS = sorfDB.find_gene('HP_1521')

    for x in sorfIDS:
        elems = sorfDB.get_sorf_infos(x)

        for x in elems:
            print(x, elems[x])
