import glob
import operator
import os
import io

import logging
from Bio import SeqIO
from collections import defaultdict, Counter

from intervaltree import IntervalTree

from assemblygraph.edge import Edge, UndirectedEdge
from assemblygraph.graph import Graph
from assemblygraph.vertex import Vertex
from database.DiamondResult import DiamondResult
from database.GeneDuplicationDB import GeneDuplicationDB

from analysis.multicombinations import strMultiCombination


from database.ModIntervalTree import ModIntervalTree
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase, MultiCombination
from utils import fileLocation



if __name__ == '__main__':

    FORMAT = '%(asctime)-15s %(user)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    d = {'user': 'mjoppich'}
    logger = logging.getLogger('homologydb')

    def log(message):
        logger.warning(message, extra=d)

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

        log("Starting Step 1")

        def removeEmptyVertices(mygraph):
            myRemoveVertexIDs = set()
            for vertexID in mygraph.vertices:

                vertex = mygraph.vertices[vertexID]

                if len(vertex.neighbors) == 0:
                    myRemoveVertexIDs.add(vertexID)

            for vertexID in myRemoveVertexIDs:
                #print("Remove: " + str(vertexID))
                mygraph.remove_vertex(vertexID)

            graph.cleanUpEmpty()

            #print(len(mygraph.vertices))
            return mygraph


        def removeBadCoverageEdges(mygraph):
            myRemoveVertexIDs = set()
            for vertexID in mygraph.vertices:

                vertex = mygraph.vertices[vertexID]

                removeEdge = list()

                for edge in vertex.neighbors:

                    edgeInfo = edge.props['info']

                    sourceSeq = edge.source.props['sequence']
                    targetSeq = edge


            for vertexID in myRemoveVertexIDs:
                mygraph.remove_vertex(vertexID)

            #print(len(mygraph.vertices))
            return mygraph


        graph = removeEmptyVertices(graph)

        """
        
        STEP 2: FIND EASY MATCHES
        
        """

        setRemoveVertexIDs = set()

        def makeIdentityScore( alignment ):
            return alignment.identity

        def makeLengthScore( alignment ):

            qseq = genomeDB.get_sequence(alignment.query.genome, alignment.query.seqid)
            sseq = genomeDB.get_sequence(alignment.subject.genome, alignment.subject.seqid)

            lengthQuery =  (len(alignment.query) / len(qseq))
            lengthSubject =(len(alignment.subject) / len(sseq))

            if lengthQuery < 0.8:
                return 0
            if lengthSubject < 0.8:
                return 0

            return (lengthQuery+lengthSubject)/2.0



        setOneVertices = set()
        for vertexID in graph.vertices:

            vertex = graph.vertices[vertexID]

            if vertex.name[1] == 'HP_0694':
                vertex.name = vertex.name

            if len(vertex.neighbors) == 1:

                targetVertex = vertex.neighbors[0].target

                # second condition is sanity check, should be always the case => unidirectional
                if len(targetVertex.neighbors) == 1 and targetVertex.neighbors[0].target == vertex:
                    diamondResult = targetVertex.neighbors[0].props['info']
                    identityScore = makeIdentityScore(diamondResult) > 0.9
                    lengthScore = makeLengthScore(diamondResult) > 0.8

                    if identityScore and lengthScore:
                        setRemoveVertexIDs.add(vertex.name)
                        setRemoveVertexIDs.add(targetVertex.name)

                        #print("Step2", vertex.name, targetVertex.name, diamondResult)
                        homolDB.addHomologyRelation(vertex.name, targetVertex.name,
                                                    {'step': "2", 'file': fileName, 'edge': (vertex.name, targetVertex.name)}
                                                    )

        for vertexID in setRemoveVertexIDs:
            graph.remove_vertex(vertexID)

        graph = removeEmptyVertices(graph)

        #print(len(graph.vertices))


        """
        
        STEP 2.1: multiple hits, but one very high scoring 
        
        """
        log("Starting Step 2.1")



        def acceptOneOfMultiple(mygraph,
                                minIdentity=0.8,
                                minQueryLength=0.85,
                                minSubjectLength=0.85,
                                edgeSortExpression=lambda x: x.props['info'].identity,
                                allowPartialLength=False,
                                allowMultiple=False,
                                betterEdgeCheck=False,
                                stepID='1ofMany'):

            sortedVerts = sorted([x for x in mygraph.vertices], key=lambda x: len(mygraph.get_vertex(x).props['sequence']),
                                 reverse=True)
            setRemoveVertexIDs = set()

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)

                if vertex.name[1] == 'HPP12_1154':
                    vertex.name=vertex.name

                for edge in sorted(vertex.neighbors, key=edgeSortExpression, reverse=True):

                    targetVertex = edge.target

                    vertexIDObj = getIDObj(edge, vertex)
                    targetVertexIDObj = getIDObj(edge, targetVertex)

                    queryLength = len(vertexIDObj) / len(vertex.props['sequence'])
                    subjectLength = len(targetVertexIDObj) / len(targetVertex.props['sequence'])

                    acceptEdge = False

                    considerEdge = False
                    if allowPartialLength:
                        considerEdge = (queryLength > minQueryLength or subjectLength > minSubjectLength) and edge.props['info'].identity > minIdentity
                        considerEdge = considerEdge and min([queryLength, subjectLength]) > 0.5
                    else:
                        considerEdge = queryLength > minQueryLength and subjectLength > minSubjectLength and edge.props['info'].identity > minIdentity

                    if considerEdge:

                        acceptEdge = True

                        for targetEdge in targetVertex.neighbors:
                            if targetEdge.props['info'].identity > edge.props['info'].identity:
                                otherVertexObj = getIDObj(targetEdge, targetEdge.target)
                                otherVertexLength = len(otherVertexObj) / len(targetEdge.target.props['sequence'])

                                if subjectLength > 0.9 and otherVertexLength > 0.9:
                                    continue

                                if otherVertexLength > subjectLength:
                                    acceptEdge = False
                                    break

                    if acceptEdge:
                        setRemoveVertexIDs.add(vertex.name)
                        setRemoveVertexIDs.add(targetVertex.name)
                        #print("acceptOneOfMultiple", minIdentity, minQueryLength, minSubjectLength, allowPartialLength, vertex.name, targetVertex.name, edge.props['info'])

                        if vertex.name[1] == 'HP_0694':
                            vertex.name=vertex.name


                        homolDB.addHomologyRelation(vertex.name, targetVertex.name,
                                                    {'step': stepID, 'file': fileName, 'edge': (vertex.name, targetVertex.name)}
                                                    )

                        if not allowMultiple:
                            break

            log("Step 2.1: Removing Vertices")
            for vertexID in setRemoveVertexIDs:
                mygraph.remove_vertex(vertexID)

            log("Step 2.1: Removing Vertices Finished")

            return mygraph


        graph = acceptOneOfMultiple(graph, allowMultiple=True)
        log("Step 2.1: Remove Empty Vertices")
        graph = removeEmptyVertices(graph)


        #print(len(graph.vertices))

        """
        
        STEP 3: One sequence, multiple sequences map
        
        """

        log("Starting Step 3")

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


        sortedVerts = sorted([x for x in graph.vertices], key=lambda x: len(graph.get_vertex(x).props['sequence']), reverse=True)
        setRemoveVertexIDs = set()

        for x in sortedVerts:
            vertex = graph.get_vertex(x)
            targetEdges = []

            tree = IntervalTree()
            baseSeq = vertex.props['sequence']
            countArray = [0] * len(baseSeq)

            allTargetsLengthMatch = True
            accumIdentity = 0.0

            if vertex.name[1] == 'jhp_0073':
                vertex.name = vertex.name


            if len(vertex.neighbors) < 2:
                continue

            for edge in sorted(vertex.neighbors, key=lambda x: x.props['info'].identity, reverse=True):

                targetVertex = edge.target
                targetSeq = targetVertex.props['sequence']

                diamondResult = edge.props['info']
                targetVertexIDObj = getIDObj(edge, targetVertex)
                vertexIDObj = getIDObj(edge, vertex)

                if targetVertexIDObj != None:
                    tree.addi( vertexIDObj.begin, vertexIDObj.end, targetVertex.name )

                    for i in range(vertexIDObj.begin-1,vertexIDObj.end):
                        countArray[ i ] += 1

                    accumIdentity += len(targetVertexIDObj) * diamondResult.identity

                    if len(targetVertexIDObj) / len(targetSeq) >= 0.7: #todo find a good value!
                        continue
                    elif len(targetVertexIDObj) / (len(targetSeq)-targetVertexIDObj.begin+1) >= 0.95: #suffix
                        continue
                    elif len(targetVertexIDObj) / (targetVertexIDObj.end) >= 0.95: #prefix
                        continue

                    allTargetsLengthMatch = False

                else:
                    raise ValueError("Problem")

            accumIdentity = accumIdentity / len(baseSeq)

            noOverlap = True
            ones = 0

            for x in countArray:
                if x > 1:
                    noOverlap = False
                    break
                elif x == 1:
                    ones += 1

            if noOverlap == False:
                continue

            possibleMatch = (ones / len(countArray) ) > 0.9 or accumIdentity > 0.5

            if vertex.name[1] == 'HP_0060':
                print(vertex)


            if possibleMatch and allTargetsLengthMatch:

                #print(vertex.name, len(baseSeq), tree)
                otherNames = set()
                for edge in vertex.neighbors:
                    otherNames.add(edge.source.name)
                    otherNames.add(edge.target.name)

                for x in otherNames:
                    setRemoveVertexIDs.add(x)

                otherNames.remove(vertex.name)

                homolDB.addCombination(vertex.name, otherNames, {'step': '3'})

        for vertexID in setRemoveVertexIDs:
            graph.remove_vertex(vertexID)

        graph = removeEmptyVertices(graph)

        #print(len(graph.vertices))


        """
        
        STEP 3.1: try to use a subset to get good coverage!
        
        SE87_04025	3	221	U063_1094	205	422	0.8540000000000001	219	31	1	7.7e-107	380.2 ('CP006888', 'U063_1094')
        SE87_04025	148	283	U063_1094	128	258	0.624	141	38	4	3.1e-39	155.6 ('CP006888', 'U063_1094') 
        SE87_04025	231	422	U063_1094	3	188	0.41100000000000003	197	100	6	1.3e-32	133.7 ('CP006888', 'U063_1094') 
        
        """

        log("Starting Step 3.1")

        def greedyBuildFromSubset(mygraph,
                                  sortingFunctionAssembly=lambda x: len(getIDObj(x, x.target)),
                                  minExplainedThreshold=0.8,
                                  allowTargetOverlaps = False
                                  ):

            sortedVerts = sorted([x for x in mygraph.vertices], key=lambda x: len(mygraph.get_vertex(x).props['sequence']), reverse=True)
            setRemoveVertexIDs = set()

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)
                targetEdges = []

                baseSeq = vertex.props['sequence']
                countArray = [0] * len(baseSeq)

                allTargetsLengthMatch = True
                accumIdentity = 0.0

                vertexTree = ModIntervalTree()
                target2tree = defaultdict(IntervalTree)
                target2vertex = dict()

                if vertex.name[1] == 'jhp_0959':
                    vertex.name = vertex.name

                if vertex.name[1] == 'HP_0091':
                    vertex.name = vertex.name

                if len(vertex.neighbors) < 2:
                    continue

                usedEdges = []

                for edge in sorted(vertex.neighbors, key=sortingFunctionAssembly, reverse=True):

                    targetVertex = edge.target
                    targetSeq = targetVertex.props['sequence']

                    diamondResult = edge.props['info']
                    targetVertexIDObj = getIDObj(edge, targetVertex)
                    vertexIDObj = getIDObj(edge, vertex)


                    part1Check = not vertexTree.overlaps(vertexIDObj.begin, vertexIDObj.end)
                    part2Check = targetVertex.name not in target2tree or not target2tree[targetVertex.name].overlaps(targetVertexIDObj.begin, targetVertexIDObj.end)

                    acceptEdge = part1Check and part2Check

                    if allowTargetOverlaps and not part1Check:
                        # mode which allows a small overlap < 15AA
                        overlaps = [len(x.intersection(vertexIDObj)) for x in vertexTree[vertexIDObj]]

                        if len(overlaps) == 0:
                            overlaps.append(0)

                        acceptOverlap = max(overlaps) < 30 and len(targetVertexIDObj) > 45
                        acceptOverlapInOther = part2Check

                        acceptEdge |= (acceptOverlap and acceptOverlapInOther)


                    if acceptEdge:

                        usedEdges.append(edge)

                        target2vertex[targetVertex.name] = targetVertex
                        vertexTree.addi(vertexIDObj.begin, vertexIDObj.end)
                        target2tree[targetVertex.name].addi(targetVertexIDObj.begin, targetVertexIDObj.end)


                if len(target2tree) == 0:
                    continue

                if vertex.name[1] in ['jhp_0054', 'jhp_0959']:
                    vertex.name = vertex.name

                vertexTree.merge_overlaps()
                totalExplained = [len(x) for x in vertexTree.all_intervals]
                explainedFraction = sum(totalExplained) / len(vertex.props['sequence'])


                if explainedFraction < 0.9 or not allowTargetOverlaps:
                    minUsedFraction = 0.9
                else:
                    minUsedFraction = 0.1

                acceptAll = True
                for targetName in target2tree:
                    used = sum([x.length() for x in target2tree[targetName]])
                    usedTargetFraction = 0.0

                    if used > 0:
                        usedTargetFraction = used / len(target2vertex[targetName].props['sequence'])

                    if used == 0 or usedTargetFraction < minUsedFraction:
                        acceptAll = False
                        break

                if len(target2tree) < 2:
                    continue

                if acceptAll == False:
                    continue


                if explainedFraction < minExplainedThreshold:
                    continue

                for x in target2vertex:
                    setRemoveVertexIDs.add(target2vertex[x].name)

                combination = set()
                for edge in usedEdges:
                    targetVertex = edge.target
                    setRemoveVertexIDs.add(targetVertex.name)
                    combination.add( targetVertex.name )

                homolDB.addCombination(vertex.name, combination, {'step': '3.1'})


            for vertexID in setRemoveVertexIDs:
                mygraph.remove_vertices_by_id(vertexID)

            mygraph = removeEmptyVertices(mygraph)

            return mygraph

        graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, 0.8)
        graph = greedyBuildFromSubset(graph, lambda x: len(getIDObj(x, x.target)), 0.8)

        #graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, minExplainedThreshold=0.5, allowTargetOverlaps=True)
        #graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, 0.55)

        #print(len(graph.vertices))

        #print(len(graph.vertices))

        """

        STEP 5: multiple sequences form a cluster

        """




        def removeDuplicateEdges( listOfEdges ):

            listOfEdgesPlus = []

            for elem in listOfEdges:
                listOfEdgesPlus.append( (elem, elem.props['info']) )

            seenResults = set()
            finalEdges = []

            for x in listOfEdgesPlus:
                if x[1] not in seenResults:
                    finalEdges.append(x[0])
                    seenResults.add(x[1])

            return finalEdges



        def makeMultiCombinations(mygraph):

            sortedVerts = sorted([x for x in mygraph.vertices],
                                 key=lambda x: len(mygraph.get_vertex(x).props['sequence']), reverse=True)

            handledVertices = set()

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)

                if vertex.name in handledVertices:
                    continue

                vertexSeq = vertex.props['sequence']

                sortingFunctionVertexEdges = lambda x: x.props['info'].identity * len(x.props['info']) #longest matches

                nextVertex = False

                if vertex.name[1] in ['HP_0733', 'HP_0732', 'jhp_0670', 'jhp_0669']:
                    vertex.name = vertex.name

                elif vertex.name[1] in ['jhp_0959', 'HPP12_1000', 'HPP12_1001', 'jhp_0958']:
                    vertex.name = vertex.name

                elif vertex.name[1] in ['jhp_0054']:
                    vertex.name = vertex.name
                else:
                    continue

                for edge in sorted(vertex.neighbors, key=sortingFunctionVertexEdges, reverse=True):

                    if nextVertex:
                        break

                    # is this a partial alignment?
                    targetVertex = edge.getOpposite(vertex)
                    targetVertexSeq = targetVertex.props['sequence']

                    vertexAlign = getIDObj(edge, vertex)
                    targetVertexAlign = getIDObj(edge, targetVertex)

                    alignedPartVertex = len(vertexAlign) / len(vertexSeq)
                    alignedPartTargetVertex = len(targetVertexAlign) / len(targetVertexSeq)

                    arePartiallyAligned = alignedPartVertex < 0.8 and alignedPartTargetVertex < 0.8 # was and - but that means that both must be partially aligned ...

                    if arePartiallyAligned:

                        allEdges = set()
                        allVertices = set()

                        combinationCreator = strMultiCombination(genomeDB)
                        combinationCreator.build_combination(edge, graph)

                        # checked all edges, not check that all elements are explained to at least 80% or whatever value
                        if not combinationCreator.valid_combination():
                            continue

                        newCombo = combinationCreator.toMultiCombination()
                        homolDB.addMultiCombination(newCombo)

                        handledVertices.union( combinationCreator.used_vertices() )

                        nextVertex = True

            for vertexID in handledVertices:
                mygraph.remove_vertex(vertexID)

            mygraph = removeEmptyVertices(mygraph)

            return mygraph


        graph = makeMultiCombinations(graph)

        """

        STEP 4: one sequence, one or multiple sequences align, accept also rather bad identity

        """

        graph = acceptOneOfMultiple(graph, 0.4, 0.8, 0.8, allowPartialLength=True, betterEdgeCheck=True,
                                    allowMultiple=False, stepID='1OfManyBad')
        graph.cleanUpEmpty()

        """
        
        STEP 6: remove hits which make no sense
        
        """
        sortedVerts = sorted([x for x in graph.vertices], key=lambda x: len(graph.get_vertex(x).props['sequence']),
                             reverse=True)
        setRemoveVertexIDs = set()
        for x in sortedVerts:
            vertex = graph.get_vertex(x)
            targetEdges = []
            baseSeq = vertex.props['sequence']
            countArray = [0] * len(baseSeq)

            for edge in sorted(vertex.neighbors, key=lambda x: x.props['info'].identity, reverse=True):

                targetVertex = edge.target
                targetSeq = targetVertex.props['sequence']

                diamondResult = edge.props['info']
                targetVertexIDObj = getIDObj(edge, targetVertex)
                vertexIDObj = getIDObj(edge, vertex)

                if targetVertexIDObj != None:
                    for i in range(vertexIDObj.begin-1,vertexIDObj.end):
                        countArray[ i ] += 1

            coverage = 0
            for x in countArray:
                if x > 0:
                    coverage += 1.0

            if vertex.name[1] == 'SE87_05730':
                pass
                #print(vertex)

            if coverage < len(baseSeq)/2:
                #not enough coverage!
                setRemoveVertexIDs.add(vertex.name)
                continue

            removeEdges = []
            for edge in vertex.neighbors:
                if edge.props['info'].identity < 0.4:
                    removeEdges.append(edge)

            for edge in removeEdges:
                idx = vertex.neighbors.index(edge)

                if idx != None and idx >= 0:
                    del vertex.neighbors[idx]

        for vertexID in setRemoveVertexIDs:
            graph.remove_vertex(vertexID)


        printGraphEdges(graph)

    homolDB.finalize()

    genomeDB.writeCSV(fileLocation+"/genome_seqs/seqs")
    homolDB.save_to_file(fileLocation + "/hpp12_hp")
