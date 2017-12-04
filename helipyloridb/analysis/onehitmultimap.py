from collections import defaultdict

from analysis.graphuser import GraphUser
from analysis.homologyresults import HomologyResult
from database.ModInterval import ModInterval
from database.ModIntervalTree import ModIntervalTree


class ManyToOneCombination(GraphUser):

    def __init__(self, graph, genomeDB, stepID='ManyToOneCombination'):

        super(ManyToOneCombination, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID

        self.used_vertex_ids = set()
        self.found_homologies = list()

    def _analyse(self):

        retRes = HomologyResult()

        sortedVerts = sorted([x for x in self.graph.vertices], key=lambda x: len(self.graph.get_vertex(x).props['sequence']), reverse=True)
        setRemoveVertexIDs = set()

        for x in sortedVerts:
            vertex = self.graph.get_vertex(x)
            targetEdges = []

            tree = ModIntervalTree()
            baseSeq = vertex.props['sequence']
            countArray = [0] * len(baseSeq)

            allTargetsLengthMatch = True
            accumIdentity = 0.0

            if vertex.name[1] in ['jhp_0073', 'HP_1551']:
                vertex.name = vertex.name

            if vertex.name in setRemoveVertexIDs:
                continue


            if len(vertex.neighbors) < 2:
                continue

            usedEdges = []
            notAccptedEdge = False

            for edge in sorted(vertex.neighbors, key=lambda x: x.props['info'].identity*x.props['info'].alignment_length, reverse=True):

                targetVertex = edge.target
                targetSeq = targetVertex.props['sequence']

                if targetVertex.name in setRemoveVertexIDs:
                    continue

                diamondResult = edge.props['info']
                targetVertexIDObj = self.getIDObj(edge, targetVertex)
                vertexIDObj = self.getIDObj(edge, vertex)

                acceptEdge = False

                if targetVertexIDObj != None:

                    if self.get_seq_fraction(edge, targetVertex) < 0.2:
                        continue

                    if len(targetVertexIDObj) / len(targetSeq) >= 0.7: #TODO either matches very long
                        acceptEdge = True
                    elif len(targetVertexIDObj) / (len(targetSeq)-targetVertexIDObj.begin+1) >= 0.95: #matches suffix
                        acceptEdge = True
                    elif len(targetVertexIDObj) / (targetVertexIDObj.end) >= 0.95: #matches prefix
                        acceptEdge = True

                    edgeInterval = ModInterval(vertexIDObj.begin, vertexIDObj.end, targetVertex.name)

                    for treeInt in tree.all_intervals:

                        interSectInt = edgeInterval.intersection(treeInt)

                        if interSectInt != None and len(interSectInt) > 15:
                            acceptEdge = False

                    if acceptEdge:

                        for i in range(vertexIDObj.begin - 1, vertexIDObj.end):
                            countArray[i] += 1

                        accumIdentity += len(targetVertexIDObj) * diamondResult.identity

                        tree.add(edgeInterval)

                        usedEdges.append( edge )

                    else:
                        # oder break?
                        notAccptedEdge = True
                        break



                else:
                    raise ValueError("Problem")

            ones = 0
            mores = 0

            for x in countArray:
                if x > 1:
                    mores += 1
                elif x == 1:
                    ones += 1

            accumIdentity = accumIdentity / len(baseSeq)
            possibleMatch = (ones / len(countArray)) > 0.9 and accumIdentity > 0.5

            if not possibleMatch:
                continue

            # if too many mores -> abort
            moreFraction = mores / len(baseSeq)

            if moreFraction > 0.1:
                continue

            # for all elements in tree, must be covered > minCoverage

            elementOverlap = defaultdict(ModIntervalTree)

            for edge in usedEdges:

                srcObj = self.getIDObj(edge, edge.source)
                elementOverlap[edge.source.name].add(srcObj)

                tgtObj = self.getIDObj(edge, edge.target)
                elementOverlap[edge.target.name].add(tgtObj)

            for nodeName in elementOverlap:
                elementOverlap[nodeName].merge_overlaps()


            invalidCoverage = False
            for nodeName in elementOverlap:
                node = self.graph.get_vertex(nodeName)
                nodeCovered = elementOverlap[nodeName].sum_intervals()
                coverage = nodeCovered / len(node.props['sequence'])

                if coverage < 0.5:
                    invalidCoverage = True
                    break


            if invalidCoverage:
                continue

            if possibleMatch and allTargetsLengthMatch:

                #print(vertex.name, len(baseSeq), tree)
                otherNames = set()
                for edge in usedEdges:
                    otherNames.add(edge.source.name)
                    otherNames.add(edge.target.name)

                for x in otherNames:
                    setRemoveVertexIDs.add(x)

                otherNames.remove(vertex.name)

                if len(otherNames) > 5:
                    print(otherNames)


                if len(otherNames) > 1:
                    retRes.combination_results.append((vertex.name, otherNames, {'step': self.stepID}))
                elif len(otherNames) == 1:

                    vName = list(otherNames)
                    retRes.homology_relations.append( (vertex.name, vName[0], {'step': self.stepID}) )
                else:
                    self.log_warn("invalid relation with 0 othernames")


        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()

        return retRes