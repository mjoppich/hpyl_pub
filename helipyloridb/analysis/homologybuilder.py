from collections import defaultdict

import math

from analysis.graphIterateAnalysis import graphCleaner
from analysis.greedymultihits import GreedyCombinationConfig, GreedyCombinationCreator
from analysis.onehithomologs import oneHitHomologs, OneHitHomologsConfig
from analysis.onehitmultimap import ManyToOneCombination
from analysis.onemultiplehits import oneMultipleHomologs, oneMultipleConfig
from analysis.removespurioushits import SpuriousEdgeRemover
from assemblygraph.graph import Graph
from assemblygraph.vertex import Vertex

from database.DiamondResult import DiamondResult
from database.GeneDuplicationDB import GeneDuplicationDB
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase
from analysis.multicombinations import MultiCombinationCreator, MultiCombinationCreatorConfig

import glob
import os


class HomologyBuilder:
    def __init__(self, basePath, inputFormat="embl", inputExtension='.gb'):

        self.basePath = basePath

        self.genomeInputExtension = inputExtension

        self.genomeDB = GenomeDB(self.basePath, fileFormat=inputFormat, fileExtension=inputExtension)
        self.homolDB = HomologyDatabase()
        self.geneDupDB = GeneDuplicationDB()

    def printResult(self, result):
        qseq = self.genomeDB.get_sequence(result.query.genome, result.query.seqid)
        sseq = self.genomeDB.get_sequence(result.subject.genome, result.subject.seqid)

        print(result.query, result.subject, result.identity, self.makeScore(result))
        print(len(qseq), qseq)
        print(len(sseq), sseq)

    def makeScore(self, result):

        iden = float(result.identity)

        qseq = self.genomeDB.get_sequence(result.query.genome, result.query.seqid)
        sseq = self.genomeDB.get_sequence(result.subject.genome, result.subject.seqid)

        length = (len(result) / len(qseq)) + (len(result) / len(sseq))

        return (4 * iden + length) / 6.0

    def getIDObj(self, edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.query

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.subject

        return None

    def getNonIDObj(self, edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.subject

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.query

        return None

    def analyse(self):

        for file in glob.glob(self.basePath + "/alignments/*.aliout"):

            query2result = defaultdict(list)
            subject2result = defaultdict(list)

            filebase = os.path.basename(file)
            afile = filebase.split('.')
            subjectGenome = afile[0]
            queryGenome = afile[1]

            fileName = filebase

            #wantedGenomes = ['AE000511', 'CP001217', 'AE001439', 'CP001173']

            wantedGenomes = None
            if wantedGenomes != None and not queryGenome in wantedGenomes:
                continue

            if wantedGenomes != None and not subjectGenome in wantedGenomes:
                continue

            if queryGenome == subjectGenome:
                continue

            self.genomeDB.loadGenome(self.basePath + "/genomes/" + queryGenome + self.genomeInputExtension)
            self.genomeDB.loadGenome(self.basePath + "/genomes/" + subjectGenome + self.genomeInputExtension)

            dupfiles = [
                self.basePath + "/alignments/" + queryGenome + "." + queryGenome + ".aliout",
                self.basePath + "/alignments/" + subjectGenome + "." + subjectGenome + ".aliout"
            ]

            canContinue = True
            for x in dupfiles:
                if not os.path.isfile(x):
                    print("Not a file", x)
                    canContinue = False

            if not canContinue:
                continue

            self.geneDupDB.load_organism(dupfiles[0],
                                         self.genomeDB)
            self.geneDupDB.load_organism(dupfiles[1],
                                         self.genomeDB)

            #print(str(self.geneDupDB))
            #print(file)

            with open(file, 'r') as infile:

                for line in infile:

                    ret = DiamondResult.from_line(line, queryGenome, subjectGenome)

                    if ret == None:
                        continue

                    if self.geneDupDB.has_gene_duplication(ret.query.genome, ret.query.seqid):
                        commonIDs = self.geneDupDB.get_gene_duplication(ret.query.genome, ret.query.seqid)
                        if ret.query.seqid != commonIDs[0]:
                            ret.query.seqid = commonIDs[0]
                            continue

                    if self.geneDupDB.has_gene_duplication(ret.subject.genome, ret.subject.seqid):
                        commonIDs = self.geneDupDB.get_gene_duplication(ret.subject.genome, ret.subject.seqid)
                        if ret.subject.seqid != commonIDs[0]:
                            ret.subject.seqid = commonIDs[0]
                            continue # todo we expect not to loose anything, but one should check that beforehand ...

                    query2result[ret.query.seqid].append(ret)
                    subject2result[ret.subject.seqid].append(ret)

            for seqid in query2result:
                allResults = query2result[seqid]
                allResults = sorted(allResults, key=lambda x: self.makeScore(x), reverse=True)
                query2result[seqid] = allResults

            for seqid in subject2result:
                allResults = subject2result[seqid]
                allResults = sorted(allResults, key=lambda x: self.makeScore(x), reverse=True)
                subject2result[seqid] = allResults


            graph = Graph()

            for seqid in query2result:
                results = query2result[seqid]

                if len(results) == 0:
                    continue

                query = results[0].query

                queryVert = Vertex(query.idtuple(), {'sequence': self.genomeDB.get_sequence(query.genome, query.seqid)})
                graph.add_vertex_if_not_exists(queryVert)

                for result in results:

                    if len(result) < 20:
                        continue

                    subjVert = Vertex(result.subject.idtuple(),
                                      {'sequence': self.genomeDB.get_sequence(result.subject.genome, result.subject.seqid)})
                    subjVert = graph.add_vertex_if_not_exists(subjVert)

                    myedge = graph.add_edge(queryVert, subjVert, {'info': result}, True)


            for vertexID in graph.vertices:
                vertex = graph.vertices[vertexID]
                vertex.neighbors = sorted(vertex.neighbors, key=lambda x: self.getNonIDObj(x, vertex).seqid)

            """

            STEP 1: REMOVE EMPTY NODES (IF EXIST)

            """
            graphClean = graphCleaner(graph, None)
            graphClean.analyse()

            """

            STEP 2: FIND EASY MATCHES

            """
            oneHitConfig = OneHitHomologsConfig()

            stepOneHits = oneHitHomologs(graph, self.genomeDB, oneHitConfig)
            homolResults = stepOneHits.analyse()
            homolResults.toDataBase(self.homolDB)

            """

            STEP 2.1: multiple hits, but one very high scoring 

            """


            """
            One of multiple. If excellent hit found, allow that one gene is homologous to many other
            
            """
            one2mulHitsConfig = oneMultipleConfig()
            one2mulHitsConfig.minQueryLength = 0.9
            one2mulHitsConfig.minSubjectLength = 0.9
            one2mulHitsConfig.allowMultiple = True

            one2mulHits = oneMultipleHomologs(graph, self.genomeDB, one2mulHitsConfig)
            homolResults = one2mulHits.analyse()
            homolResults.toDataBase(self.homolDB)


            """
            
            One of multiple, allow non-excellent hits
            
            """
            one2mulHitsConfig = oneMultipleConfig()
            one2mulHitsConfig.allowMultiple = False

            one2mulHits = oneMultipleHomologs(graph, self.genomeDB, one2mulHitsConfig)
            homolResults = one2mulHits.analyse()
            homolResults.toDataBase(self.homolDB)


            """

            STEP 3: One sequence, multiple sequences map

            """
            many2one = ManyToOneCombination(graph, self.genomeDB)
            retRes = many2one.analyse()
            retRes.toDataBase(self.homolDB)

            """

            STEP 3.1: try to use a subset to get good coverage!

            """

            greedyConfig = GreedyCombinationConfig()
            greedyConfig.sortingFunctionAssembly = lambda x: x.props['info'].identity
            greedyCreator = GreedyCombinationCreator(graph, self.genomeDB, greedyConfig)
            retRes = greedyCreator.analyse()
            retRes.toDataBase(self.homolDB)



            greedyConfig = GreedyCombinationConfig()
            greedyConfig.minExplainedThreshold=0.5
            greedyConfig.allowTargetOverlaps=True

            greedyCreator = GreedyCombinationCreator(graph, self.genomeDB, greedyConfig)
            retRes = greedyCreator.analyse()
            retRes.toDataBase(self.homolDB)

            greedyConfig = GreedyCombinationConfig()
            greedyConfig.minExplainedThreshold=0.55

            greedyCreator = GreedyCombinationCreator(graph, self.genomeDB, greedyConfig)
            retRes = greedyCreator.analyse()
            retRes.toDataBase(self.homolDB)

            """

            STEP 5: multiple sequences form a cluster

            """

            mulCombAnalysisConfig = MultiCombinationCreatorConfig()
            mulCombAnalysis = MultiCombinationCreator(graph, self.genomeDB, mulCombAnalysisConfig)
            mulCombResult = mulCombAnalysis.analyse()
            mulCombResult.toDataBase(self.homolDB)

            """

            STEP 4: one sequence, one or multiple sequences align, accept also rather bad identity

            """

            omConfig = oneMultipleConfig()
            omConfig.allowPartialLength = True
            omConfig.betterEdgeCheck = True
            omConfig.allowMultiple = False
            omConfig.minIdentity = 0.4
            omConfig.minQueryLength = 0.8
            omConfig.minSubjectLength = 0.8

            def checkEdge(config, edge, source, target):

                queryLength = config.get_seq_fraction(edge, source)
                subjectLength = config.get_seq_fraction(edge, target)

                considerEdge = queryLength > config.minQueryLength
                considerEdge = considerEdge and subjectLength > config.minSubjectLength
                considerEdge = considerEdge and edge.props['info'].identity > config.minIdentity

                considerEdge = considerEdge and min([queryLength, subjectLength]) > 0.5

                return considerEdge

            omConfig.considerEdgeFunc = checkEdge

            omAnalysis = oneMultipleHomologs(graph, self.genomeDB, omConfig)
            retRes = omAnalysis.analyse()
            retRes.toDataBase(self.homolDB)

            """

            extremely long sequences > 500!

            """


            def checkEdgeLong(config, edge, source, target):

                edgeInfo = edge.props['info']

                minSeqLength = min([len(source.props['sequence']), len(target.props['sequence'])])

                queryLength = config.get_seq_fraction(edge, source)
                subjectLength = config.get_seq_fraction(edge, target)

                if edgeInfo.identity * minSeqLength > 500:

                    if edgeInfo.evalue < math.pow(10, -90):

                        if queryLength > config.minQueryLength and subjectLength > config.minQueryLength:
                            return True

                return False

            omConfig.minQueryLength = 0.6
            omConfig.minSubjectLength = 0.6
            omConfig.considerEdgeFunc = checkEdgeLong

            omAnalysis = oneMultipleHomologs(graph, self.genomeDB, omConfig)
            retRes = omAnalysis.analyse()
            retRes.toDataBase(self.homolDB)

            """

            STEP 6: remove hits which make no sense

            """

            edgeRemover = SpuriousEdgeRemover(graph, self.genomeDB)
            edgeRemover.analyse()


            """
            
            Some relations may have been hidden by combinations
            
            """
            oneHitConfig = OneHitHomologsConfig(minIDScore=0.8, minLengthScore=0.7)

            stepOneHits = oneHitHomologs(graph, self.genomeDB, oneHitConfig)
            homolResults = stepOneHits.analyse()
            homolResults.toDataBase(self.homolDB)


            """

            STEP 7: Mention leftovers

            """

            def printEdge(edge):
                print(edge.source.name, edge.target.name, edge.props['info'])

            def printGraphEdges(mygraph):

                sortedVerts = sorted([x for x in mygraph.vertices],
                                     key=lambda x: len(mygraph.get_vertex(x).props['sequence']),
                                     reverse=True)

                seenDiamondInfos = set()
                for x in sortedVerts:

                    vertex = mygraph.get_vertex(x)

                    for edge in vertex.neighbors:

                        diamondInfo = edge.props.get('info', None)

                        if diamondInfo == None:
                            continue

                        if diamondInfo not in seenDiamondInfos:
                            printEdge(edge)
                            seenDiamondInfos.add(diamondInfo)

            printGraphEdges(graph)



        allDupRelations = self.geneDupDB.get_gene_relations()
        allDupRelations.toDataBase(self.homolDB)

        self.homolDB.finalize()

        self.genomeDB.writeCSV(self.basePath + "/genome_seqs/seqs")

        return self.homolDB