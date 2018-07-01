import argparse
import re
import shlex
import sys, os
from io import StringIO

from database.PfamResultDB import PfamResultDB

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../helipyloridb")

from flask import Flask, jsonify, request, redirect, url_for, send_from_directory
import json
import pprint
from collections import defaultdict
import requests

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from database.OperonDB import OperonDB
from database.SORFDB import SORFDB
from database.TSSDB import TSSDB
from database.XRefDatabase import XRefDatabase
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from utils import fileLocation

from Bio.Align.Applications import ClustalOmegaCommandline

from Bio import SeqIO, AlignIO
import subprocess
from flask_cors import CORS

dataurl = str(os.path.dirname(os.path.realpath(__file__))) + "/../../" + 'frontend/src/static/'

dataurl = os.path.abspath(dataurl)

print("Loading website from", dataurl)

app = Flask(__name__, static_folder=dataurl, static_url_path='/static')
CORS(app)

app.config['DEBUG'] = False
app.config['UPLOAD_FOLDER'] = ""


# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']


# genomDB.writeBLASTfastas(fileLocation + "/genomes")
# genomDB.writePFAMfastas('/mnt/c/Users/mjopp/Desktop/genomdb')
@app.route('/')
def root():

    retFile = 'index.html'

    return app.send_static_file(retFile)

@app.route('/test', methods=['GET', 'POST'])
def test():
    return "<html><body>heliPyloriDB Server v0.01</body></html>", 200, None


@app.route('/help', methods=['GET', 'POST'])
def help():
    res = "<html><body><ul>"

    for x in [rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static']:
        res += "<li>" + str(x) + "</li>"

    res += "</body></html>"

    return res, 200, None


@app.route('/findHomCluster/<geneid>')
def findHomCluster(geneid):
    if geneid == None or not geneid in homDB.homologies:
        return app.make_response((jsonify({'error': 'invalid homID'}), 400, None))

    homIDs = homDB.findHomologyForGeneID(geneid)

    jsonResult = {'homs': homIDs}

    return app.make_response((jsonify(jsonResult), 200, None))


@app.route('/homcluster/<homID>')
def getHomCluster(homID):
    if homID == None or not homID in homDB.homologies:
        return app.make_response((jsonify({'error': 'invalid homID'}), 400, None))

    homCluster = homDB.get_homology_cluster(homID)

    jsonResult = {'elements': homCluster}

    return app.make_response((jsonify(jsonResult), 200, None))


@app.route('/homcluster', methods=['GET', 'POST'])
def homcluster():
    homID = request.values.get('homid', None)

    return getHomCluster(homID)


@app.route('/alignment', methods=['GET', 'POST'])
def getAlignment():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    if not 'homid' in alignReq:
        return app.make_response((jsonify({'error': 'must include homid'}), 400, None))

    homID = alignReq['homid']
    alignOrgs = alignReq['organisms'] if 'organisms' in alignReq else homDB.get_all_organisms()

    return returnAlignments([homID], alignOrgs)


@app.route('/clustalign', methods=['GET', 'POST'])
def getHomClusterAndAlignment():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    if not 'genes' in alignReq and not 'homs' in alignReq:
        return app.make_response((jsonify({'error': 'must include gene or homid'}), 400, None))

    gene2homid = {}
    if 'genes' in alignReq:
        for geneID in alignReq['genes']:
            gene2homid[geneID] = homDB.findHomologyForGeneID(geneID)

    if 'homs' in alignReq:
        for homid in alignReq['homs']:
            gene2homid[homid] = [homid]

    alignOrgs = alignReq['organisms'] if (
                'organisms' in alignReq and len(alignReq['organisms']) > 0) else homDB.get_all_organisms()

    return returnAlignments(gene2homid, alignOrgs)


@app.route('/alignments', methods=['GET', 'POST'])
def getAlignments():
    alignReq = request.get_json(force=True, silent=True)

    if alignReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    if not 'homids' in alignReq:
        return app.make_response((jsonify({'error': 'must include homid'}), 400, None))

    homID = alignReq['homids']
    alignOrgs = alignReq['organisms'] if (
                'organisms' in alignReq and len(alignReq['organisms']) > 0) else homDB.get_all_organisms()

    return returnAlignments({homID: homID}, alignOrgs)


def returnAlignments(homIDs, alignOrgs):
    print("Preparing Alignment for ", homIDs)
    print("Preparing Alignment for ", alignOrgs)

    jsonResult = defaultdict(list)

    for refID in homIDs:

        for homID in homIDs[refID]:
            homCluster = homDB.get_cluster(homID)

            if homCluster == None:
                continue

            alignSeqs = []

            for org in alignOrgs:
                if org in homCluster:
                    for seqid in homCluster[org]:
                        alignSeqs.append((org, seqid))

            seqRecords = []
            seqID2Element = {}
            for org, seqid in alignSeqs:
                genSeq = genomDB.get_sequence(org, seqid)

                seqRecID = "_".join([org, seqid])

                seqID2Element[seqRecID] = genomDB.get_element(org, seqid)
                seq = SeqRecord(Seq(genSeq, generic_dna), id=seqRecID, description="")

                seqRecords.append(seq)

            if len(seqRecords) > 0:
                outfasta = open(tmpfolder + '/out.fasta', 'w')
                SeqIO.write(seqRecords, outfasta, "fasta")
                outfasta.close()

                inmsa = tmpfolder + "/in.msa"

                clustalomega_cline = ClustalOmegaCommandline(cmd=clustalobin, infile=outfasta.name, outfile=inmsa, force=True,
                                                             outfmt='fa', verbose=True, auto=True)
                print(clustalomega_cline)
                output = subprocess.getoutput([str(clustalomega_cline)])
                print("Clustalomega finished")

                alignment = []
                with open(inmsa, 'r') as fin:
                    alignment = AlignIO.read(fin, "fasta")

                clusterMSA = []

                tssList = set()
                operonsList = set()
                sorfList = set()
                allPfamRes = []

                for seqr in alignment:

                    ida = seqr.id.split('_', 1)

                    if ida[0] == 'AE000511':

                        if opDB.find_gene(ida[1], None) != None:

                            inOperons = opDB.find_gene(ida[1])

                            for operon in inOperons:
                                operonsList.add(operon)

                        if tssDB.find_gene(ida[1], None) != None:

                            inTSS = tssDB.find_gene(ida[1])

                            for tssid in inTSS:
                                tssList.add(tssid)

                        if sorfDB.find_gene(ida[1], None) != None:
                            inSORF = sorfDB.find_gene(ida[1])

                            for sorfid in inSORF:
                                sorfList.add(sorfid)

                    orgID = ida[0]
                    genID = ida[1]

                    pfamResIDs = pfamDB.find_gene(orgID, genID)
                    pfamRes = pfamDB.get_pfam_infos(pfamResIDs)

                    allPfamRes += pfamRes

                    genomeEntry = seqID2Element[seqr.id].toJSON().copy()
                    foundXRefs = xrefDB.make_infos(ida[1])

                    genomeEntry['alignment'] = str(seqr.seq)
                    genomeEntry['alignmentNT'] = makeNTCoAlign(seqr.seq, genomeEntry)
                    genomeEntry['xrefs'] = foundXRefs

                    clusterMSA.append(
                        genomeEntry)  # {'org': ida[0], 'seqid': ida[1], 'alignment': str(seqr.seq), 'xrefs': foundXRefs} )

                operonsInfo = []
                for operon in operonsList:
                    operonsInfo.append(opDB.get_operon_infos(operon))

                tssInfo = []
                for tssid in tssList:
                    tssInfo.append(tssDB.get_tss_infos(tssid))

                sorfInfo = []
                for sorfid in sorfList:
                    sorfInfo.append(sorfDB.get_sorf_infos(sorfid))

                for x in allPfamRes:
                    print(x)

                jsonResult[refID].append(
                    {'msa': clusterMSA, 'homid': homID, 'TSS': tssInfo, 'OPERONS': operonsInfo, 'SORFS': sorfInfo,
                     'PFAMS': allPfamRes})

    print(jsonResult)

    return app.make_response((jsonify(jsonResult), 200, None))


def makeNTCoAlign(alignAA, genEntry):
    ntseq = genEntry['seqNT']

    usedCodons = []
    for istart in range(0, len(ntseq), 3):
        usedCodons.append(ntseq[istart:istart + 3])

    alignmentNT = ""
    alignmentAA = ""
    seqPos = 0

    try:

        for seqChar in str(alignAA):
            if seqChar == '-':
                alignmentNT += "---"
                alignmentAA += "---"
            else:
                alignmentNT += usedCodons[seqPos]
                alignmentAA += seqChar + "--"
                seqPos += 1

    except:
        print("Error in retrieving NT align")
        print(genEntry)
        alignmentAA = None

    return (alignmentNT, alignmentAA)


@app.route('/autocomplete', methods=['GET', 'POST'])
def findID():
    jsonResult = {}
    jsonResult['locus_tag'] = []
    jsonResult['HOMS'] = []

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 3:
        return app.make_response((jsonify(jsonResult), 200, None))

    reMatch = re.compile(searchWord)

    for homid in homDB.get_all_ids():
        if reMatch.match(homid):
            jsonResult['HOMS'].append(homid)

            if len(jsonResult['HOMS']) > 100:
                break

    for orgname in homDB.get_all_organisms():

        for protname in genomDB.get_sequences_for_genome(orgname):

            if reMatch.match(protname):
                jsonResult['locus_tag'].append(protname)

                if len(jsonResult['locus_tag']) > 100:
                    break

    allResults = []

    for x in jsonResult:
        for elem in jsonResult[x]:
            allResults.append({'name': elem, 'type': x})

    return app.make_response((jsonify(allResults), 200, None))


@app.route('/orgautocomplete', methods=['GET', 'POST'])
def orgac():
    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 3:
        return app.make_response((jsonify({}), 200, None))

    reMatch = re.compile(searchWord)

    allorgs = homDB.get_all_organisms()
    allorgs = [{'name': x, 'id': x} for x in allorgs if reMatch.match(x)]  # beautify names

    return app.make_response((jsonify(allorgs), 200, None))


@app.route('/defaultorgs', methods=['GET', 'POST'])
def defaultorgs():
    allorgs = homDB.get_all_organisms()
    allorgs = [{'name': x, 'id': x} for x in allorgs]  # beautify names

    allorgs = [{'name': x, 'id': x} for x in ['AE000511', 'CP001217', 'AE001439']]  # beautify names

    return app.make_response((jsonify(allorgs), 200, None))


@app.route('/organisms', methods=['GET', 'POST'])
def get_organisms():
    allorgs = homDB.get_all_organisms()
    allorgs = [{'name': x, 'id': x} for x in allorgs]  # beautify names

    return app.make_response((jsonify(allorgs), 200, None))


@app.route('/stats', methods=['GET', 'POST'])
def make_simple_stats():
    allorgs = homDB.get_all_organisms()

    jsonRes = {
        'org_count': len(allorgs),
        'hom_count': len(homDB.get_hom_ids()),
        'comb_count': len(homDB.get_comb_ids()),
        'mul_comb_count': len(homDB.get_mulcombs())
    }

    return app.make_response((jsonify(jsonRes), 200, None))


@app.route('/swissmodel/query', methods=['POST'])
def make_swissmodel_query():
    searchWords = request.get_json(force=True, silent=True)

    print("SWISSMODEL " + str(searchWords))

    if searchWords == None or searchWords.get('uniprot', None) == None:
        return app.make_response((jsonify({"error": "no uniprot id given", 'req': searchWords}), 400, None))

    uniprotIDs = searchWords.get('uniprot')

    retRes = {}

    for uniprotID in uniprotIDs:
        try:
            r = requests.get(
                'https://swissmodel.expasy.org/repository/uniprot/' + uniprotID + '.json?provider=swissmodel')

            if r.status_code == 200:
                retRes[uniprotID] = r.json()

        except:
            pass

    return app.make_response((jsonify(retRes), 200, None))


# homDB = HomologyDatabase.loadFromFile(fileLocation + "/hpp12_hp")

homDB = None
xrefDB = None
opDB = None
tssDB = None
sorfDB = None
pfamDB = None
genomDB = None
tmpfolder = "/tmp/"
clustalobin = None


def start_app_from_args(args):
    global homDB
    global genomDB
    global xrefDB
    global opDB
    global sorfDB
    global pfamDB
    global tssDB
    global tmpfolder
    global clustalobin


    tmpfolder = args.tmp
    clustalobin = args.clustalo.name

    homDB = HomologyDatabase.loadFromFile(args.databases + "/homdb/" + "/hpdb_full_new")
    xrefDB = XRefDatabase(args.databases + "/homdb/" + "/hpdb_full_xref")
    opDB = OperonDB.from_cs_operons(args.databases + "/sharma/operons.xlsx")
    tssDB = TSSDB.from_cs_tss(args.databases + "/sharma/tss.xlsx")
    sorfDB = SORFDB.from_cs_sorfs(args.databases + "/sharma/sorfs.xlsx")
    pfamDB = PfamResultDB.from_folder(args.databases + "/pfam/")

    genomDB = GenomeDB(args.genomes, loadAll=False)

    for orgname in homDB.get_all_organisms():
        genomDB.loadGenome(orgname)


def getCLParser():
    parser = argparse.ArgumentParser(description='Start hpylDB Data Server', add_help=False)
    parser.add_argument('-g', '--genomes', type=str, required=True)
    parser.add_argument('-d', '--databases', type=str, required=True)
    parser.add_argument('--tmp', type=str, required=True)
    parser.add_argument('--clustalo', type=argparse.FileType('r'), required=False, default="/usr/bin/clustalo")
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    return parser


if __name__ == '__main__':

    parser = getCLParser()

    args = parser.parse_args()

    for x in args.__dict__:
        print(x, args.__dict__[x])

    start_app_from_args(args)

    print("Starting Flask on port", args.port)
    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)


def gunicorn_start(genomes, databases, tmpdir):
    parser = getCLParser()

    argstr = "--genomes {genomesdir} --databases {databases} --tmp {tmpdir}".format(genomesdir=genomes, databases=databases, tmpdir=tmpdir)

    print("Starting app with")
    print(argstr)

    args = parser.parse_args(shlex.split(argstr))

    start_app_from_args(args)

    return app