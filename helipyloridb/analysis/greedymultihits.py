from collections import defaultdict

from analysis.graphuser import GraphUser
from analysis.homologyresults import HomologyResult
from database.ModIntervalTree import ModIntervalTree


class GreedyCombinationConfig:

    def __init__(self):

        self.minExplainedThreshold = 0.8
        self.allowTargetOverlaps = False

        self.sortingFunctionAssembly = None

class GreedyCombinationCreator(GraphUser):

    def __init__(self, graph, genomeDB, config, stepID='GreedyCombinationCreator'):

        super(GreedyCombinationCreator, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID
        self.config = config

        self.used_vertex_ids = set()
        self.found_homologies = list()


    def _analyse(self):

        returnResults = HomologyResult()

        sortedVerts = sorted([x for x in self.graph.vertices], key=lambda x: len(self.graph.get_vertex(x).props['sequence']),
                             reverse=True)
        setRemoveVertexIDs = set()

        for x in sortedVerts:
            vertex = self.graph.get_vertex(x)
            targetEdges = []

            baseSeq = vertex.props['sequence']
            countArray = [0] * len(baseSeq)

            allTargetsLengthMatch = True
            accumIdentity = 0.0

            vertexTree = ModIntervalTree()
            target2tree = defaultdict(ModIntervalTree)
            target2vertex = dict()

            if vertex.name[1] == 'jhp_0959':
                vertex.name = vertex.name

            if vertex.name[1] == 'HP_0091':
                vertex.name = vertex.name

            if len(vertex.neighbors) < 2:
                continue

            usedEdges = []

            if self.config.sortingFunctionAssembly == None:
                sortingFunctionAssembly = lambda x: len(self.getIDObj(x, x.target))
            else:
                sortingFunctionAssembly = self.config.sortingFunctionAssembly


            for edge in sorted(vertex.neighbors, key=sortingFunctionAssembly, reverse=True):

                targetVertex = edge.target
                targetSeq = targetVertex.props['sequence']

                diamondResult = edge.props['info']
                targetVertexIDObj = self.getIDObj(edge, targetVertex)
                vertexIDObj = self.getIDObj(edge, vertex)

                part1Check = not vertexTree.overlaps(vertexIDObj.begin, vertexIDObj.end)
                part2Check = targetVertex.name not in target2tree or not target2tree[targetVertex.name].overlaps(
                    targetVertexIDObj.begin, targetVertexIDObj.end)

                acceptEdge = part1Check and part2Check

                if self.config.allowTargetOverlaps and not part1Check:
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

            if explainedFraction < 0.9 or not self.config.allowTargetOverlaps:
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

            if explainedFraction < self.config.minExplainedThreshold:
                continue

            for x in target2vertex:
                setRemoveVertexIDs.add(target2vertex[x].name)

            combination = set()
            for edge in usedEdges:
                targetVertex = edge.target
                setRemoveVertexIDs.add(targetVertex.name)
                combination.add(targetVertex.name)

            returnResults.combination_results.append((vertex.name, combination, {'step': '3.1'}))

        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()


        return returnResults