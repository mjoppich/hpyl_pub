import glob
import re
from collections import defaultdict

import os
from openpyxl import load_workbook


class PfamResultDB:

    def __init__(self):

        self.pfams = defaultdict(lambda: defaultdict(list))
        self.org2pfamresid = defaultdict(lambda: defaultdict(list))



    def find_gene(self, organism, geneID, defaultval=list()):

        retIDs = []

        for pfamresid in self.org2pfamresid.get(organism, []):

            pfamres = self.org2pfamresid[organism][pfamresid]

            if geneID == pfamres['GENE_ID']:
                retIDs.append(pfamresid)

        return defaultval if len(retIDs) == 0 else retIDs

    def get_pfam_info(self, pfamresid):

        if not pfamresid in self.pfams:
            return None

        pfamData = self.pfams[pfamresid]

        return pfamData

    def get_pfam_infos(self, pfamresids):

        res = []
        for pfamresid in pfamresids:
            res.append(self.get_pfam_info(pfamresid))

        return res


    @classmethod
    def from_folder(cls, fileslocation="/mnt/c/ownCloud/data/hpyloriDB/pfam/res/"):

        pfamDB = PfamResultDB()


        curID = 0

        for file in glob.glob(fileslocation + "*.pfam"):

            fileBasename = os.path.basename(file)
            orgID = os.path.splitext(os.path.splitext(fileBasename)[0])[0]

            with open(file, 'r') as fin:

                for line in fin:

                    if len(line) == 0 or line.startswith('#'):
                        continue

                    line = line.strip()
                    aline = line.split()

                    if len(aline) == 0:
                        continue

                    hasActiveSites = len(aline) == 16

                    pfamResID = "PFAMRES" + str(curID)
                    curID += 1

                    # <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_residues>

                    seqID = aline[0]
                    alignStart = int(aline[1])
                    alignEnd = int(aline[2])
                    envStart = int(aline[3])
                    envEnd = int(aline[4])

                    pfamID = aline[5]
                    pfamName = aline[6]

                    pfamType = aline[7]

                    pfamHMMStart = int(aline[8])
                    pfamHMMEnd = int(aline[9])
                    pfamHMMLength = int(aline[10])

                    pfamBitScore = float(aline[11])
                    pfamEvalue = float(aline[12])
                    pfamSignificant = int(aline[13])

                    pfamClan = aline[14]

                    if pfamClan == 'No_clan':
                        pfamClan = None

                    pfamActiveSites = []

                    if hasActiveSites:

                        res = re.findall(r'\[(.*?)\]', aline[15])

                        actSites = set()

                        for actSite in res:
                            ints = actSite.split(",")
                            ints = [int(x) for x in ints]

                            for x in ints:
                                actSites.add(x)

                        pfamActiveSites = actSites


                    pfamDB.org2pfamresid[orgID][pfamResID] = {'GENE_ID': seqID}

                    pfamRes= {
                        'GENE_ID': seqID,
                        'GENE_START': alignStart,
                        'GENE_END': alignEnd,
                        'ENV_START': envStart,
                        'ENV_END': envEnd,
                        'PFAMID': pfamID,
                        'PFAMNAME': pfamName,
                        'PFAM_TYPE': pfamType,
                        'HMM_START': pfamHMMStart,
                        'HMM_END': pfamHMMEnd,
                        'HMM_LENGTH': pfamHMMLength,
                        'PFAM_BIT_SCORE': pfamBitScore,
                        'PFAM_EVALUE': pfamEvalue,
                        'PFAM_SIGNIFICANT': pfamSignificant
                    }

                    if pfamClan:
                        pfamRes['PFAM_CLAN'] = pfamClan

                    if hasActiveSites:
                        pfamRes['PFAM_ACTIVE_SITES'] = pfamActiveSites

                    pfamDB.pfams[pfamResID] = pfamRes

        return pfamDB


if __name__ == '__main__':

    pfamDB = PfamResultDB.from_folder()

    pfamResIDs = pfamDB.find_gene('AE000511','HP_0626')

    for x in pfamResIDs:
        print(pfamDB.get_pfam_info(x))
