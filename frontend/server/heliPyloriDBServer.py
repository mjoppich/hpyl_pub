import re
import sys
from io import StringIO

from flask import Flask, jsonify, request, redirect, url_for, send_from_directory
import os
import json
import pprint
from collections import defaultdict

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from database.OperonDB import OperonDB
from database.TSSDB import TSSDB
from database.XRefDatabase import XRefDatabase
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation

from Bio.Align.Applications import ClustalOmegaCommandline

from Bio import SeqIO, AlignIO
import subprocess
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

app.config['DEBUG'] = False
app.config['UPLOAD_FOLDER'] = ""

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp12_hp")
genomDB = GenomeDB(fileLocation + "/genomes", loadAll=False)
xrefDB = XRefDatabase()
opDB = OperonDB.from_cs_operons()
tssDB = TSSDB.from_cs_tss()


for orgname in homDB.get_all_organisms():
    genomDB.loadGenome(orgname)


@app.route('/test', methods=['GET', 'POST'])
def test():

    return "<html><body>heliPyloriDB Server v0.01</body></html>", 200, None


@app.route('/help', methods=['GET', 'POST'])
def help():
    res = "<html><body><ul>"

    for x in [rule.rule for rule in app.url_map.iter_rules() if rule.endpoint !='static']:
        res += "<li>"+str(x)+"</li>"

    res +="</body></html>"

    return res, 200, None


@app.route('/findHomCluster/<geneid>')
def findHomCluster(geneid):

    if geneid == None or not geneid in homDB.homologies:
        return app.make_response((jsonify( {'error': 'invalid homID'} ), 400, None))

    homIDs = homDB.findHomologyForGeneID(geneid)

    jsonResult = {'homs': homIDs}

    return app.make_response((jsonify( jsonResult ), 200, None))


@app.route('/homcluster/<homID>')
def getHomCluster(homID):

    if homID == None or not homID in homDB.homologies:
        return app.make_response((jsonify( {'error': 'invalid homID'} ), 400, None))


    homCluster = homDB.get_homology_cluster(homID)

    jsonResult = {'elements': homCluster}

    return app.make_response((jsonify( jsonResult ), 200, None))


@app.route('/homcluster', methods=['GET', 'POST'])
def homcluster():

    homID = request.values.get('homid', None)

    return getHomCluster(homID)

@app.route('/alignment', methods=['GET', 'POST'])
def getAlignment():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))


    if not 'homid' in alignReq:
        return app.make_response((jsonify( {'error': 'must include homid'} ), 400, None))

    homID = alignReq['homid']
    alignOrgs = alignReq['organisms'] if 'organisms' in alignReq else homDB.get_all_organisms()

    return returnAlignments([homID], alignOrgs)


@app.route('/clustalign', methods=['GET', 'POST'])
def getHomClusterAndAlignment():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    if not 'genes' in alignReq:
        return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


    gene2homid = {}
    for geneID in alignReq['genes']:
        gene2homid[geneID] = homDB.findHomologyForGeneID(geneID)


    alignOrgs = alignReq['organisms'] if 'organisms' in alignReq else homDB.get_all_organisms()

    return returnAlignments(gene2homid, alignOrgs)

@app.route('/alignments', methods=['GET', 'POST'])
def getAlignments():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'homids' in alignReq:
        return app.make_response((jsonify( {'error': 'must include homid'} ), 400, None))

    homID = alignReq['homids']
    alignOrgs = alignReq['organisms'] if 'organisms' in alignReq else homDB.get_all_organisms()

    return returnAlignments({homID: homID}, alignOrgs)


def returnAlignments(homIDs, alignOrgs):

    jsonResult = defaultdict(list)

    for refID in homIDs:

        for homID in homIDs[refID]:
            homCluster = homDB.get_cluster(homID)

            alignSeqs = []

            for org in alignOrgs:
                if org in homCluster:

                    for seqid in homCluster[org]:

                        alignSeqs.append( (org, seqid) )


            seqRecords = []
            seqID2Element = {}
            for org, seqid in alignSeqs:

                genSeq = genomDB.get_sequence(org, seqid)

                seqRecID = "_".join([org, seqid])

                seqID2Element[seqRecID] = genomDB.get_element(org, seqid)
                seq = SeqRecord(Seq(genSeq, generic_dna), id=seqRecID, description="")

                seqRecords.append(seq)

            outfasta = open('./tmps/out.fasta', 'w')
            SeqIO.write(seqRecords, outfasta, "fasta")
            outfasta.close()

            inmsa = "./tmps/in.msa"

            clustalomega_cline = ClustalOmegaCommandline(infile=outfasta.name, outfile=inmsa,force=True, outfmt='st', verbose=True, auto=True)
            print(clustalomega_cline)

            output = subprocess.getoutput( [str(clustalomega_cline)] )

            alignment = []
            with open(inmsa, 'r') as fin:
                alignment = AlignIO.read(fin, "stockholm")


            clusterMSA = []

            tssList = set()
            operonsList = set()

            for seqr in alignment:

                ida = seqr.id.split('_', 1)

                if ida[0] == 'AE000511':

                    if opDB.find_gene(ida[1], None) != None:

                        inOperons = opDB.find_gene(ida[1])

                        for operon in inOperons:
                            operonsList.add( operon )

                    if tssDB.find_gene(ida[1], None) != None:

                        inTSS = tssDB.find_gene(ida[1])

                        for tssid in inTSS:
                            tssList.add(tssid)


                genomeEntry = seqID2Element[seqr.id].toJSON().copy()
                foundXRefs = xrefDB.make_infos(ida[1])

                genomeEntry['alignment'] = str(seqr.seq)

                genomeEntry['alignmentNT'] = makeNTCoAlign(seqr.seq, genomeEntry)
                genomeEntry['xrefs'] = foundXRefs

                clusterMSA.append( genomeEntry)#{'org': ida[0], 'seqid': ida[1], 'alignment': str(seqr.seq), 'xrefs': foundXRefs} )

            operonsInfo = []
            for operon in operonsList:
                operonsInfo.append(opDB.get_operon_infos(operon))

            tssInfo = []
            for tssid in tssList:
                tssInfo.append(tssDB.get_tss_infos(tssid))

            jsonResult[refID].append({'msa': clusterMSA, 'homid': homID, 'TSS': tssInfo, 'OPERONS': operonsInfo})

    return app.make_response((jsonify( jsonResult ), 200, None))

def makeNTCoAlign(alignAA, genEntry):

    ntseq = genEntry['seqNT']

    usedCodons = []
    for istart in range(0, len(ntseq), 3):
        usedCodons.append(ntseq[istart:istart+3])

    alignmentNT = ""
    alignmentAA = ""
    seqPos = 0

    for seqChar in str(alignAA):
        if seqChar == '-':
            alignmentNT += "---"
            alignmentAA += "---"
        else:
            alignmentNT += usedCodons[seqPos]
            alignmentAA += seqChar + "--"
            seqPos += 1

    return (alignmentNT, alignmentAA)



@app.route('/autocomplete', methods=['GET', 'POST'])
def findID():

    jsonResult = {}
    jsonResult['genomes'] = []
    jsonResult['proteins'] = []

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 3:
        return app.make_response((jsonify( jsonResult ), 200, None))



    reMatch = re.compile(searchWord)

    for orgname in homDB.get_all_organisms():

        if reMatch.match(orgname):
            jsonResult['genomes'].append(orgname)

    for orgname in homDB.get_all_organisms():

        for protname in genomDB.get_sequences_for_genome(orgname):

            if reMatch.match(protname):
                jsonResult['proteins'].append(protname)

    return app.make_response((jsonify( jsonResult ), 200, None))



if __name__ == '__main__':

   print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint !='static'])

   app.run(threaded=True)
