import glob
import operator
import os
import io

import logging
from Bio import SeqIO
from collections import defaultdict, Counter

from intervaltree import IntervalTree

from analysis.graphIterateAnalysis import graphCleaner
from analysis.greedymultihits import GreedyCombinationConfig, GreedyCombinationCreator
from analysis.onehithomologs import oneHitHomologs
from analysis.onehitmultimap import ManyToOneCombination
from analysis.onemultiplehits import oneMultipleHomologs, oneMultipleConfig
from analysis.removespurioushits import SpuriousEdgeRemover
from assemblygraph.edge import Edge, UndirectedEdge
from assemblygraph.graph import Graph
from assemblygraph.vertex import Vertex
from database.DiamondResult import DiamondResult
from database.GeneDuplicationDB import GeneDuplicationDB

from analysis.multicombinations import MultiCombinationCreator, MultiCombinationCreatorConfig


from database.ModIntervalTree import ModIntervalTree
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase, MultiCombination
from utils import fileLocation


if __name__ == '__main__':



    genomeDB = GenomeDB(fileLocation)
    homolDB = HomologyDatabase()
    geneDupDB = GeneDuplicationDB()

    def printResult(result):
        qseq = genomeDB.get_sequence(result.query.genome, result.query.seqid)
        sseq = genomeDB.get_sequence(result.subject.genome, result.subject.seqid)

        print(result.query, result.subject, result.identity, makeScore(result))
        print(len(qseq), qseq)
        print(len(sseq), sseq)

    def makeScore(result):

        iden = float(result.identity)

        qseq = genomeDB.get_sequence(result.query.genome, result.query.seqid)
        sseq = genomeDB.get_sequence(result.subject.genome, result.subject.seqid)

        length = (len(result) / len(qseq)) + (len(result) / len(sseq))

        return (4*iden + length) / 6.0


    def getIDObj(edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.query

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.subject

        return None


    def getNonIDObj(edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.subject

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.query

        return None


    for file in glob.glob(fileLocation + "/alignments/*.aliout"):

        query2result = defaultdict(list)
        subject2result = defaultdict(list)

        filebase = os.path.basename(file)
        afile = filebase.split('.')
        subjectGenome = afile[0]
        queryGenome = afile[1]

        fileName = filebase


        wantedGenomes = ['AE000511', 'CP001217', 'AE001439']

        if not queryGenome in wantedGenomes:
            continue

        if not subjectGenome in wantedGenomes:
            continue

        if queryGenome == subjectGenome:
            continue

        genomeDB.loadGenome(fileLocation + "/genomes/" + queryGenome + ".gb")
        genomeDB.loadGenome(fileLocation + "/genomes/" + subjectGenome + ".gb")

        geneDupDB.load_organism(fileLocation + "/alignments/" + queryGenome + "." + queryGenome + ".aliout", genomeDB)
        geneDupDB.load_organism(fileLocation + "/alignments/" + subjectGenome + "." + subjectGenome + ".aliout", genomeDB)

        print(str(geneDupDB))

        print(file)

        with open(file, 'r') as infile:

            for line in infile:

                ret = DiamondResult.from_line(line, queryGenome, subjectGenome)

                if ret == None:
                    continue

                if geneDupDB.has_gene_duplication(ret.query.genome, ret.query.seqid):
                    commonIDs = geneDupDB.get_gene_duplication(ret.query.genome, ret.query.seqid)
                    ret.query.seqid = commonIDs[0]

                if geneDupDB.has_gene_duplication(ret.subject.genome, ret.subject.seqid):
                    commonIDs = geneDupDB.get_gene_duplication(ret.subject.genome, ret.subject.seqid)
                    ret.subject.seqid = commonIDs[0]


                query2result[ret.query.seqid].append(ret)
                subject2result[ret.subject.seqid].append(ret)

        for seqid in query2result:
            allResults = query2result[seqid]
            allResults = sorted(allResults, key=lambda x: makeScore(x), reverse=True)
            query2result[seqid] = allResults

        for seqid in subject2result:
            allResults = subject2result[seqid]
            allResults = sorted(allResults, key=lambda x: makeScore(x), reverse=True)
            subject2result[seqid] = allResults

        usedIDs = set()

        graph = Graph()

        for seqid in query2result:
            results = query2result[seqid]

            if len(results) == 0:
                continue

            query = results[0].query

            queryVert = Vertex(query.idtuple(), {'sequence': genomeDB.get_sequence(query.genome, query.seqid)})
            graph.add_vertex_if_not_exists(queryVert)

            for result in results:

                if len(result) < 20:
                    continue

                subjVert = Vertex(result.subject.idtuple(), {'sequence': genomeDB.get_sequence(result.subject.genome, result.subject.seqid)})
                subjVert = graph.add_vertex_if_not_exists(subjVert)

                myedge = graph.add_edge(queryVert, subjVert, {'info': result}, True)



        #print(len(graph.vertices))

        for vertexID in graph.vertices:
            vertex = graph.vertices[vertexID]
            vertex.neighbors = sorted(vertex.neighbors, key=lambda x: getNonIDObj(x, vertex).seqid)

        """
        
        STEP 1: REMOVE EMPTY NODES (IF EXIST)
        
        """
        graphClean = graphCleaner(graph, None)
        graphClean.analyse()


        """
        
        STEP 2: FIND EASY MATCHES
        
        """
        stepOneHits = oneHitHomologs(graph, genomeDB)
        homolResults = stepOneHits.analyse
        homolResults.toDataBase(homolDB)


        """
        
        STEP 2.1: multiple hits, but one very high scoring 
        
        """

        one2mulHitsConfig = oneMultipleConfig()

        one2mulHits = oneMultipleHomologs(graph, genomeDB, one2mulHitsConfig)
        homolResults = one2mulHits.analyse()
        homolResults.toDataBase(homolDB)

        """
        
        STEP 3: One sequence, multiple sequences map
        
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



        def printGraph(mygraph):

            sortedVerts = sorted([x for x in mygraph.vertices], key=lambda x: len(mygraph.get_vertex(x).props['sequence']),
                                 reverse=True)

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)
                targetEdges = []

                #print(vertex)

                for x in vertex.neighbors:

                    tv = mygraph.get_vertex(x.target)
                    add = ""
                    seqname = "NONE!!!"
                    if tv != None:
                        add = tv.props['sequence']
                        seqname = tv.name

                    #print(x.props['info'], seqname, add)


        many2one = ManyToOneCombination(graph, genomeDB)
        retRes = many2one.analyse()
        retRes.toDataBase(homolDB)


        """
        
        STEP 3.1: try to use a subset to get good coverage!
                
        """

        greedyConfig = GreedyCombinationConfig()
        greedyConfig.sortingFunctionAssembly = lambda x: x.props['info'].identity

        greedyCreator = GreedyCombinationCreator(graph, genomeDB, greedyConfig)
        retRes = greedyCreator.analyse()
        retRes.toDataBase(homolDB)


        greedyConfig = GreedyCombinationConfig()
        greedyCreator = GreedyCombinationCreator(graph, genomeDB, greedyConfig)
        retRes = greedyCreator.analyse()
        retRes.toDataBase(homolDB)

        #graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, minExplainedThreshold=0.5, allowTargetOverlaps=True)
        #graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, 0.55)

        #print(len(graph.vertices))

        #print(len(graph.vertices))

        """

        STEP 5: multiple sequences form a cluster

        """

        mulCombAnalysisConfig = MultiCombinationCreatorConfig()
        mulCombAnalysis = MultiCombinationCreator(graph, genomeDB, mulCombAnalysisConfig)
        mulCombResult = mulCombAnalysis.analyse()
        mulCombResult.toDataBase(homolDB)


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

        omAnalysis = oneMultipleHomologs(graph, genomeDB, omConfig)
        retRes = omAnalysis.analyse()

        retRes.toDataBase(homolDB)

        """
        
        STEP 6: remove hits which make no sense
        
        """

        edgeRemover = SpuriousEdgeRemover(graph, genomeDB)
        edgeRemover.analyse()


        printGraphEdges(graph)

    homolDB.finalize()

    genomeDB.writeCSV(fileLocation+"/genome_seqs/seqs")
    homolDB.save_to_file(fileLocation + "/hpp12_hp")
