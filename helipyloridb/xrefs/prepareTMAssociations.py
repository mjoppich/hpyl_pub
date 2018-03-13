from porestat.utils.DataFrame import DataFrame, DataRow

from utils import fileLocation
from xrefs.Synfile import Synfile

from collections import defaultdict

syndict = {}
for line in open(fileLocation + "/tm/output/synfile.map"):

    aline = line.strip().split(": ")
    syndict[aline[1]] = Synfile(aline[0])


outDataFrame = DataFrame()
outDataFrame.addColumns(['XREF', 'GOID'])

for tmFile in ["/tm/output/PfamFamily.index", "/tm/output/interpro.index"]:

    with open(fileLocation + tmFile, 'r') as fin:

        pfamID2GO = defaultdict(set)

        for line in fin:

            aline = line.split("\t")

            pfamID = aline[0].split(".")[0]

            synpos = aline[1]
            asynpos = synpos.split(':')

            synFileID = asynpos[0]
            synRowID = int(asynpos[1])

            synFile = syndict[synFileID]
            synID = synFile.line2syn[synRowID]

            foundSyn = synFile.mSyns[synID]

            goID = foundSyn.id.replace('_', ':')

            pfamID2GO[pfamID].add(goID)

        for pfamID in pfamID2GO:
            for goid in pfamID2GO[pfamID]:

                drow = DataRow.fromDict({
                    'XREF': pfamID,
                    'GOID': goid
                })

                outDataFrame.addRow(drow)

outDataFrame.export(fileLocation + "hpp12_hp_xref_add")
