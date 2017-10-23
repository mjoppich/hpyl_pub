import glob
import os
import io
from Bio import SeqIO
from collections import defaultdict, Counter

from intervaltree import IntervalTree

from assemblygraph.edge import Edge
from assemblygraph.graph import Graph
from assemblygraph.vertex import Vertex
from database.DiamondResult import DiamondResult
from database.ModInterval import ModInterval
from database.ModIntervalTree import ModIntervalTree
from database.genomedb import GenomeDB
from database.homologydb import HomologyDatabase, MultiCombination
from utils import fileLocation



if __name__ == '__main__':

    genomeDB = GenomeDB(fileLocation)
    homolDB = HomologyDatabase()

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


    for file in glob.glob(fileLocation + "/genomes/alis/*.aliout"):

        query2result = defaultdict(list)
        subject2result = defaultdict(list)

        filebase = os.path.basename(file)
        afile = filebase.split('.')
        subjectGenome = afile[0]
        queryGenome = afile[1]


        wantedGenomes = ['AE000511', 'CP001217', 'AE001439']

        if not queryGenome in wantedGenomes:
            continue

        if not subjectGenome in wantedGenomes:
            continue

        genomeDB.loadGenome(fileLocation + "/genomes/" + queryGenome + ".gb")
        genomeDB.loadGenome(fileLocation + "/genomes/" + subjectGenome + ".gb")

        with open(file, 'r') as infile:

            for line in infile:

                ret = DiamondResult.from_line(line, queryGenome, subjectGenome)

                if ret == None:
                    continue

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
                subjVert = Vertex(result.subject.idtuple(), {'sequence': genomeDB.get_sequence(result.subject.genome, result.subject.seqid)})
                subjVert = graph.add_vertex_if_not_exists(subjVert)

                graph.add_edge(queryVert, subjVert, {'info': result}, True)

        #print(len(graph.vertices))

        """
        
        STEP 1: REMOVE EMPTY NODES (IF EXIST)
        
        """

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

            lengthQuery =  (len(alignment) / len(qseq))
            lengthSubject =(len(alignment) / len(sseq))

            if lengthQuery < 0.8:
                return 0
            if lengthSubject < 0.8:
                return 0

            return (lengthQuery+lengthSubject)/2.0



        setOneVertices = set()
        for vertexID in graph.vertices:

            vertex = graph.vertices[vertexID]

            if len(vertex.neighbors) == 1:

                targetVertex = vertex.neighbors[0].target

                # second condition is sanity check, should be always the case => unidirectional
                if len(targetVertex.neighbors) == 1 and targetVertex.neighbors[0].target == vertex:
                    setRemoveVertexIDs.add(vertex.name)
                    setRemoveVertexIDs.add(targetVertex.name)


                    diamondResult = targetVertex.neighbors[0].props['info']
                    identityScore = makeIdentityScore(diamondResult) > 0.9
                    lengthScore = makeLengthScore(diamondResult) > 0.8

                    if identityScore and lengthScore:
                        #print("Step2", vertex.name, targetVertex.name, diamondResult)
                        homolDB.addHomologyRelation(vertex.name, targetVertex.name, {'step': "2"})

        for vertexID in setRemoveVertexIDs:
            graph.remove_vertex(vertexID)

        graph = removeEmptyVertices(graph)

        #print(len(graph.vertices))


        """
        
        STEP 2.1: multiple hits, but one very high scoring 
        
        """

        def getIDObj(edge, vertex):

            diamondResult = edge.props['info']

            if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
                return diamondResult.query

            if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
                return diamondResult.subject

            return None


        def acceptOneOfMultiple(mygraph, minIdentity=0.8, minQueryLength=0.85, minSubjectLength=0.85, allowPartialLength=False, allowMultiple=False,betterEdgeCheck=False, stepID='1ofMany'):

            sortedVerts = sorted([x for x in mygraph.vertices], key=lambda x: len(mygraph.get_vertex(x).props['sequence']),
                                 reverse=True)
            setRemoveVertexIDs = set()

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)

                if vertex.name[1] == 'SE87_00145':
                    pass
                    #print(vertex.name)

                for edge in sorted(vertex.neighbors, key=lambda x: x.props['info'].identity, reverse=True):

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

                        if vertex.name[1] == 'HPP12_1404':
                            #print(vertex)
                            pass

                        homolDB.addHomologyRelation(vertex.name, targetVertex.name, {'step': stepID})

                        if not allowMultiple:
                            break

            for vertexID in setRemoveVertexIDs:
                mygraph.remove_vertex(vertexID)

            return mygraph


        graph = acceptOneOfMultiple(graph, allowMultiple=True)
        graph = removeEmptyVertices(graph)


        #print(len(graph.vertices))

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

            if vertex.name[1] == 'HP_0424':
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

        def greedyBuildFromSubset(mygraph, sortingFunctionAssembly=lambda x: len(getIDObj(x, x.target)), minExplainedThreshold=0.8):

            sortedVerts = sorted([x for x in mygraph.vertices], key=lambda x: len(mygraph.get_vertex(x).props['sequence']), reverse=True)
            setRemoveVertexIDs = set()

            for x in sortedVerts:
                vertex = mygraph.get_vertex(x)
                targetEdges = []

                baseSeq = vertex.props['sequence']
                countArray = [0] * len(baseSeq)

                allTargetsLengthMatch = True
                accumIdentity = 0.0

                vertexTree = IntervalTree()
                target2tree = defaultdict(IntervalTree)
                target2vertex = dict()

                if vertex.name[1] == 'HP_0764':
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

                    if not vertexTree.overlaps(vertexIDObj.begin, vertexIDObj.end) and not target2tree[targetVertex.name].overlaps(targetVertexIDObj.begin, targetVertexIDObj.end):

                        usedEdges.append(edge)

                        target2vertex[targetVertex.name] = targetVertex
                        vertexTree.addi(vertexIDObj.begin, vertexIDObj.end)
                        target2tree[targetVertex.name].addi(targetVertexIDObj.begin, targetVertexIDObj.end)


                if len(target2tree) == 0:
                    continue

                acceptAll = True
                totalExplained = 0
                for targetName in target2tree:
                    used = sum([x.length() for x in target2tree[targetName]])
                    totalExplained += used

                    if used / len(target2vertex[targetName].props['sequence']) < 0.8:
                        acceptAll = False
                        break

                if len(target2tree) < 2:
                    continue

                if acceptAll == False:
                    continue

                if totalExplained/len(vertex.props['sequence']) < minExplainedThreshold:
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
                mygraph.remove_vertex(vertexID)

            mygraph = removeEmptyVertices(mygraph)

            return mygraph

        graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, 0.8)
        graph = greedyBuildFromSubset(graph, lambda x: len(getIDObj(x, x.target)), 0.8)
        graph = greedyBuildFromSubset(graph, lambda x: x.props['info'].identity, 0.55)

        #print(len(graph.vertices))

        """
        
        STEP 4: one sequence, one or multiple sequences align, accept also rather bad identity
        
        """

        graph = acceptOneOfMultiple(graph, 0.4, 0.8, 0.8, allowPartialLength=True, betterEdgeCheck=True, allowMultiple=False, stepID='1OfManyBad')
        graph.cleanUpEmpty()

        #print(len(graph.vertices))

        """

        STEP 5: multiple sequences form a cluster

        """


        def checkIntervalOverlap(setOfOverlaps, baseInterval):

            if setOfOverlaps == None:
                return True

            sourceOverlapOk = True
            if not len(setOfOverlaps) == 0:
                for x in setOfOverlaps:
                    interSect = x.intersection(baseInterval)
                    if interSect == None:
                        continue

                    if len(interSect) > 10:
                        sourceOverlapOk = False
                        break

            return sourceOverlapOk


        def getAllEdges(vertex, knownEdges=set(), seenVertices=set()):

            if vertex == None:
                return set()

            if len(vertex.neighbors) == 0:
                return set()

            foundEdges = knownEdges

            for edge in vertex.neighbors:
                foundEdges.add(edge)

            seenVertices.add(vertex)

            for edge in vertex.neighbors:

                if vertex == edge.source and edge.target not in seenVertices:
                    seenVertices.add(edge.target)

                    (newEdges, newVertices) = getAllEdges(edge.target, foundEdges, seenVertices)
                    foundEdges = foundEdges.union(newEdges)
                    seenVertices = seenVertices.union(newVertices)

                elif edge.source not in seenVertices:

                    seenVertices.add(edge.source)
                    (newEdges, newVertices) = getAllEdges(edge.source, foundEdges, seenVertices)
                    foundEdges = foundEdges.union(newEdges)
                    seenVertices = seenVertices.union(newVertices)

            return (foundEdges, seenVertices)

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

                sortingFunctionAssembly = lambda x: x.props['info'].identity

                nextVertex = False

                if vertex.name[1] in ['HP_0733', 'HP_0732', 'jhp_0670', 'jhp_0669']:
                    vertex.name = vertex.name

                for edge in sorted(vertex.neighbors, key=sortingFunctionAssembly, reverse=True):

                    if nextVertex:
                        break

                    # is this a partial alignment?
                    targetVertex = edge.getOpposite(vertex)
                    targetVertexSeq = targetVertex.props['sequence']

                    vertexAlign = getIDObj(edge, vertex)
                    targetVertexAlign = getIDObj(edge, targetVertex)

                    alignedPartVertex = len(vertexAlign) / len(vertexSeq)
                    alignedPartTargetVertex = len(targetVertexAlign) / len(targetVertexSeq)

                    bothPartiallyAligned = alignedPartVertex < 0.8 and alignedPartTargetVertex < 0.8

                    if bothPartiallyAligned:

                        allEdges = set()
                        allVertices = set()

                        # find zusammenhangskomponente
                        (allEdges, allVertices) = getAllEdges(vertex, allEdges, allVertices)
                        (allEdges2, allVertices2) = getAllEdges(targetVertex, allEdges, allVertices)

                        allEdges = allEdges.union(allEdges2)
                        allEdges = removeDuplicateEdges(allEdges)

                        usedRelations = set()
                        coveredRegions = defaultdict(ModIntervalTree)

                        def testKeyFunc(x):
                            return x.props['info'].identity

                        for foundEdge in sorted(allEdges, key=testKeyFunc, reverse=True):

                            feSource = foundEdge.source
                            feTarget = foundEdge.target

                            feSourceAlign = getIDObj(foundEdge, feSource)
                            feTargetAlign = getIDObj(foundEdge, feTarget)

                            # is there any overlap in source ?
                            overlapsSource = coveredRegions[feSourceAlign.idtuple()][feSourceAlign]
                            overlapsTarget = coveredRegions[feTargetAlign.idtuple()][feTargetAlign]

                            overlapsSourceOK = checkIntervalOverlap(overlapsSource, feSourceAlign)
                            overlapsTargetOK = checkIntervalOverlap(overlapsTarget, feTargetAlign)

                            if not overlapsSourceOK or not overlapsTargetOK:
                                continue  # cannot use this edge

                            usedRelations.add(
                                (
                                    ModInterval(feSourceAlign.begin, feSourceAlign.end, feSourceAlign.idtuple()),
                                    ModInterval(feTargetAlign.begin, feTargetAlign.end, feTargetAlign.idtuple())
                                )
                            )

                            coveredRegions[feSourceAlign.idtuple()].add(feSourceAlign)
                            coveredRegions[feTargetAlign.idtuple()].add(feTargetAlign)

                            # merge any overlapping intervals
                            coveredRegions[feSourceAlign.idtuple()].merge_overlaps(newinttype=ModInterval)
                            coveredRegions[feTargetAlign.idtuple()].merge_overlaps(newinttype=ModInterval)

                        # checked all edges, not check that all elements are explained to at least 80% or whatever value
                        allWellExplained = True

                        for seqid in coveredRegions:

                            seqTree = coveredRegions[seqid]
                            seqTree.merge_overlaps()

                            seqTreeCovered = sum([len(x) for x in seqTree])

                            sequence = genomeDB.get_sequence(seqid[0], seqid[1])

                            explainedFraction = seqTreeCovered / len(sequence)

                            if explainedFraction < 0.8:
                                allWellExplained = False
                                break

                        if allWellExplained:

                            originalRelations = usedRelations
                            usedRelations = list(usedRelations)

                            while len(usedRelations) > 0:

                                currentCombo = set()
                                currentCombo.add(usedRelations[0])

                                lastCurrentComboSize = len(currentCombo) - 1

                                while lastCurrentComboSize != len(currentCombo):

                                    lastCurrentComboSize = len(currentCombo)

                                    usedVertices = set()
                                    for rel in currentCombo:
                                        usedVertices.add(rel[0].data)
                                        usedVertices.add(rel[1].data)

                                    for rel in usedRelations:
                                        if rel[0].data in usedVertices:
                                            currentCombo.add(rel)
                                        if rel[1].data in usedVertices:
                                            currentCombo.add(rel)

                                for x in currentCombo:
                                    usedRelations.remove(x)

                                newCombo = MultiCombination()

                                for relation in currentCombo:

                                    handledVertices.add(relation[0].data)
                                    handledVertices.add(relation[1].data)

                                    newCombo.addMatch(relation[0], relation[1])

                                homolDB.addMultiCombination(newCombo)

                            nextVertex = True

            mygraph = removeEmptyVertices(mygraph)

            return mygraph


        graph = makeMultiCombinations(graph)

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
    homolDB.printCombinations()

    print("HomolDB")

    genomeDB.writeCSV(fileLocation+"/genome_seqs/seqs")

    with open(fileLocation+"/hpp12_hp", 'w') as outfile:
        outfile.write(str(homolDB))


    exit(0)


