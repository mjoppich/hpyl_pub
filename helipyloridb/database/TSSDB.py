import re
from collections import defaultdict

from openpyxl import load_workbook


class TSSDB:

    def __init__(self):

        self.tss = defaultdict(list)
        self.tss_props = defaultdict(lambda: dict())


    def find_gene(self, geneID, defaultval=list()):

        retOperonIDs = []

        for x in self.tss:

            if geneID == self.tss[x][0]:
                retOperonIDs.append(x)

        return defaultval if len(retOperonIDs) == 0 else retOperonIDs

    def get_tss_infos(self, tssID):

        if not tssID in self.tss:
            return None

        tssData = self.tss[tssID]
        tssProps = self.tss_props[tssID]


        #locusTag, tss, strand, tssSeq

        tssBase= {
            'TSSID': tssID,
            'LOCUSTAG': tssData[0],
            'TSS': tssData[1],
            'STRAND': tssData[2],
            'SEQ': tssData[3],
            'PROPS': tssProps
        }

        return tssBase


    @classmethod
    def from_cs_tss(cls, filelocation="/mnt/c/ownCloud/data/hpyloriDB/sharma/tss.xlsx"):

        tssDB = TSSDB()

        wb = load_workbook(filelocation)
        ws = wb.active

        curTSSID = 0

        primaryCol = 4
        secondaryCol = 5
        internalCol = 6
        antisenseCol = 7
        enrichedCol = 8
        locationCol = 9
        putative_sRNACol = 10
        putative_asRNACol = 11
        seqCol = 12


        for row in ws:
            if row[0].row < 9:
                continue

            locusTag = row[2].value
            locusTag = locusTag[0:2]+'_' + re.sub(r"[^0-9]", "", locusTag[2:])

            strand = row[1].value
            tss = row[0].value

            tssSeq = row[seqCol].value


            tssProps = []

            if row[primaryCol].value == 1:
                tssProps.append('primary')

            if row[secondaryCol].value == 1:
                tssProps.append('secondary')

            if row[internalCol].value == 1:
                tssProps.append('internal')

            if row[antisenseCol].value == 1:
                tssProps.append('antisense')

            if row[enrichedCol].value == 1:
                tssProps.append('enriched')

            if row[locationCol].value == 1:
                tssProps.append('location')

            if row[putative_sRNACol].value == 1:
                tssProps.append('putative sRNA')

            if row[putative_asRNACol].value == 1:
                tssProps.append('putative asRNA')


            curTss = 'TSS' + str(curTSSID)
            curTSSID += 1

            tssDB.tss[curTss] = (locusTag, tss, strand, tssSeq)
            tssDB.tss_props[curTss] = tssProps

            #print(locusTag, tss, strand, tssSeq, tssProps)
        return tssDB


if __name__ == '__main__':

    tssDB = TSSDB.from_cs_tss()

    tssIDs = tssDB.find_gene('HP_1523')

    for x in tssIDs:
        print(tssDB.get_tss_infos(x))
