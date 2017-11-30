from analysis.graphuser import GraphUser
from analysis.homologyresults import HomologyResult
from database.ModIntervalTree import ModIntervalTree


class ManyToOneCombination(GraphUser):

    def __init__(self, graph, genomeDB, stepID='ManyToOneCombination'):

        super(ManyToOneCombination, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID

        self.used_vertex_ids = set()
        self.found_homologies = list()

    def analyse(self):

        retRes = HomologyResult()

        sortedVerts = sorted([x for x in self.graph.vertices], key=lambda x: len(graph.get_vertex(x).props['sequence']), reverse=True)
        setRemoveVertexIDs = set()

        for x in sortedVerts:
            vertex = self.graph.get_vertex(x)
            targetEdges = []

            tree = ModIntervalTree()
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
                targetVertexIDObj = self.getIDObj(edge, targetVertex)
                vertexIDObj = self.getIDObj(edge, vertex)

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

                retRes.combination_results.append((vertex.name, otherNames, {'step': '3'}))

        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()

        return retRes