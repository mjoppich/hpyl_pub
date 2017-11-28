import operator
from collections import defaultdict
from intervaltree import IntervalTree

from assemblygraph.graph import Graph
from database.ModInterval import ModInterval
from database.ModIntervalTree import ModIntervalTree
from database.homologydb import MultiCombination


class GraphUser:

    def getIDObj(self, edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.query

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.subject

        return None

class strMultiCombination(GraphUser):

    def __init__(self, genomeDB):

        super(strMultiCombination, self).__init__()

        self.genomeDB = genomeDB

        self.possible_edges = None
        self.considered_edges = None
        self.used_relations = None
        self.covered_regions = None
        self.used_edges = None

        self.reset()

        self.ok_explained_fraction = 0.7

        self.well_max_gap_fraction = 0.1
        self.well_explained_fraction = 0.8
        self.well_explained_max_gap = 200
        self.well_explained_covered_abs = 125

        self.valid_combination_fraction = 0.6

    def checkIntervalOverlap(self, setOfOverlaps, baseInterval):

        overlaps = self.interval_overlap(setOfOverlaps, baseInterval)
        sampledOverlaps = [x for x in overlaps if len(x) < 10]

        return len(overlaps) == len(sampledOverlaps)

    def toMultiCombination(self):

        retComb = MultiCombination()

        for edge in self.used_edges:

            edgeInfo = edge.props['info']

            retComb.addMatch( edgeInfo.query, edgeInfo.subject )

        return retComb

    def used_vertices(self):

        vertexIDs = set()

        for edge in self.used_edges:
            vertexIDs.add(edge.source.name)
            vertexIDs.add(edge.target.name)

        return vertexIDs


    def interval_overlap(self, setOfOverlaps, baseInterval):

        overlaps = []

        if setOfOverlaps == None:
            return overlaps

        if not len(setOfOverlaps) == 0:
            for x in setOfOverlaps:
                interSect = x.intersection(baseInterval)

                if interSect == None:
                    continue

                overlaps.append(interSect)

        return overlaps


    def getAllEdges(self, vertex, knownEdges=set(), seenVertices=set(), minOverlapLength=20):

        if vertex == None:
            return set()

        if len(vertex.neighbors) == 0:
            return set()

        foundEdges = knownEdges

        for edge in vertex.neighbors:
            foundEdges.add(edge)

        seenVertices.add(vertex)

        for edge in vertex.neighbors:

            if len(edge.props['info']) < minOverlapLength:
                continue

            if vertex == edge.source and edge.target not in seenVertices:
                seenVertices.add(edge.target)

                (newEdges, newVertices) = self.getAllEdges(edge.target, foundEdges, seenVertices, minOverlapLength)
                foundEdges = foundEdges.union(newEdges)
                seenVertices = seenVertices.union(newVertices)

            elif edge.source not in seenVertices:

                seenVertices.add(edge.source)
                (newEdges, newVertices) = self.getAllEdges(edge.source, foundEdges, seenVertices, minOverlapLength)
                foundEdges = foundEdges.union(newEdges)
                seenVertices = seenVertices.union(newVertices)

        return (foundEdges, seenVertices)


    def reset(self):
        self.possible_edges = None
        self.considered_edges = set()
        self.used_relations = set()
        self.used_edges = set()

        self.covered_regions = defaultdict(ModIntervalTree)



    def build_combination(self, startEdge, graph):
        # add closest edges gradually ...

        self.reset()

        combinationGraph = Graph()

        nextCandidateEdge = startEdge

        while nextCandidateEdge != None:

            feSource = nextCandidateEdge.source
            feTarget = nextCandidateEdge.target

            feSourceAlign = self.getIDObj(nextCandidateEdge, feSource)
            feTargetAlign = self.getIDObj(nextCandidateEdge, feTarget)

            overlapsSource = self.covered_regions[feSourceAlign.idtuple()][feSourceAlign]
            overlapsTarget = self.covered_regions[feTargetAlign.idtuple()][feTargetAlign]

            overlapsSourceOK = self.checkIntervalOverlap(overlapsSource, feSourceAlign)
            overlapsTargetOK = self.checkIntervalOverlap(overlapsTarget, feTargetAlign)

            if not overlapsSourceOK or not overlapsTargetOK:
                useEdge = False  # cannot use this edge
            else:
                useEdge = True

            if useEdge:
                combinationGraph.add_vertex_if_not_exists(nextCandidateEdge.source)
                combinationGraph.add_vertex_if_not_exists(nextCandidateEdge.target)

                combinationGraph.add_edge(nextCandidateEdge.source, nextCandidateEdge.target,
                                          {'info': nextCandidateEdge.props['info'], 'origEdge': nextCandidateEdge}, undirected=True)

                self.used_relations.add(
                    (
                        ModInterval(feSourceAlign.begin, feSourceAlign.end, [nextCandidateEdge.props['info']]),
                        ModInterval(feTargetAlign.begin, feTargetAlign.end, [nextCandidateEdge.props['info']])
                    )
                )

                self.covered_regions[feSourceAlign.idtuple()].add(feSourceAlign)
                self.covered_regions[feTargetAlign.idtuple()].add(feTargetAlign)

                def mergeLists(x,y):
                    if x == None:
                        x = []
                    elif y == None:
                        y = []

                    return x+y

                # merge any overlapping intervals
                self.covered_regions[feSourceAlign.idtuple()].merge_overlaps(newinttype=ModInterval, data_initializer=list(), data_reducer=mergeLists)
                self.covered_regions[feTargetAlign.idtuple()].merge_overlaps(newinttype=ModInterval, data_initializer=list(), data_reducer=mergeLists)

                self.used_edges.add(nextCandidateEdge)

            self.considered_edges.add(nextCandidateEdge)

            well_explained = self.well_explained()

            (nextCandidateEdge, explainedBases) = self.calculateNextCandidateEdge(graph, combinationGraph)

            if well_explained and explainedBases < 0:
                break

            if not well_explained and explainedBases < -5.0:
                break

    def calculateNextCandidateEdge(self, originalGraph, combinationGraph):

        # from combination graph get all nodes
        allNodes = combinationGraph.get_vertices()

        considerableEdges = dict()

        # from originalGraph get all edges which have not yet been added to considered_edges
        for edge in originalGraph.get_edges(
                decision_function=lambda edge: edge.source.name in allNodes or edge.target.name in allNodes):

            # must not be used yet
            if edge in self.considered_edges:
                continue

            # must share node
            if not edge.source.name in allNodes and not edge.target.name in allNodes:
                continue

            # calculate how much an edge would improve the coverage of elements in combinationGraph
            improvedBases = 0

            edgeInfo = edge.props['info']
            fIdentity = edgeInfo.identity

            srcSeq = edge.source.props['sequence']
            tgtSeq = edge.target.props['sequence']

            srcSeqFract = len(self.getIDObj(edge, edge.source)) / len(srcSeq)
            tgtSeqFract = len(self.getIDObj(edge, edge.target)) / len(tgtSeq)

            minCovFract = 0.5
            edgeCheck = srcSeqFract > minCovFract or tgtSeqFract > minCovFract

            if not edgeCheck:
                continue

            if edgeInfo.subject.idtuple() in allNodes:
                overlaps = self.covered_regions[edgeInfo.subject.idtuple()][edgeInfo.subject]
                intersects = self.interval_overlap(overlaps, edgeInfo.subject)
                overlap = sum([len(x) for x in intersects])

                improvement = len(edgeInfo.subject) * fIdentity - overlap
                improvedBases += improvement

            if edgeInfo.query.idtuple() in allNodes:
                overlaps = self.covered_regions[edgeInfo.query.idtuple()][edgeInfo.query]
                intersects = self.interval_overlap(overlaps, edgeInfo.query)
                overlap = sum([len(x) for x in intersects])

                improvement = len(edgeInfo.query) * fIdentity - overlap
                improvedBases += improvement

            considerableEdges[edge] = improvedBases

        # return best edge
        sortedEdges = sorted(considerableEdges.items(), key=operator.itemgetter(1), reverse=True)

        if len(sortedEdges) == 0:
            return (None, -1)

        return sortedEdges[0]

    def valid_combination(self):

        for seqid in self.covered_regions:

            seqTree = self.covered_regions[seqid]
            seqTree.merge_overlaps()

            seqTreeCovered = sum([len(x) for x in seqTree])
            sequence = self.genomeDB.get_sequence(seqid[0], seqid[1])

            sequenceInterval = ModInterval(1, len(sequence))
            sequenceTree = ModIntervalTree([sequenceInterval])

            for x in seqTree:
                sequenceTree.chop(x.begin, x.end)

            uncoveredGapIntervals = sorted([x for x in sequenceTree.all_intervals], key=lambda x: x.begin)
            uncoveredGaps = [len(x) for x in uncoveredGapIntervals]
            maxGap = max(uncoveredGaps) if len(uncoveredGaps) > 0 else -1

            explainedFraction = seqTreeCovered / len(sequence)

            if explainedFraction < self.valid_combination_fraction:
                return False

        return True

    def well_explained(self):
        allWellExplained = True

        for seqid in self.covered_regions:

            seqTree = self.covered_regions[seqid]
            seqTree.merge_overlaps()

            seqTreeCovered = sum([len(x) for x in seqTree])
            sequence = self.genomeDB.get_sequence(seqid[0], seqid[1])

            sequenceInterval = ModInterval(1, len(sequence))
            sequenceTree = ModIntervalTree([sequenceInterval])

            for x in seqTree:
                sequenceTree.chop(x.begin, x.end)

            uncoveredGapIntervals = sorted([x for x in sequenceTree.all_intervals], key=lambda x: x.begin)
            uncoveredGaps = [len(x) for x in uncoveredGapIntervals]
            maxGap = max(uncoveredGaps) if len(uncoveredGaps) > 0 else -1

            explainedFraction = seqTreeCovered / len(sequence)

            if explainedFraction > self.well_explained_fraction:
                continue

            if self.ok_explained_fraction < explainedFraction and explainedFraction < self.well_explained_fraction:

                gapCond = maxGap > min(self.well_explained_max_gap, self.well_max_gap_fraction*len(sequence))

                # borderCond: either no gap or if gap, first or last interval
                borderCond = maxGap < 0 or (maxGap > 0 and (uncoveredGaps.index(maxGap) == 0 or uncoveredGaps.index(maxGap) == len(uncoveredGaps)-1))

                if borderCond:
                    gapCond = maxGap > min(self.well_explained_max_gap, 2 * self.well_max_gap_fraction * len(sequence))

                if gapCond:
                    allWellExplained = False
                    return False

            if explainedFraction < self.ok_explained_fraction:
                allWellExplained = False
                return False

        return allWellExplained
