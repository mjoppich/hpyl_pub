from analysis.graphuser import GraphUser


class SpuriousEdgeRemover(GraphUser):

    def __init__(self, graph, genomeDB, stepID='SpuriousEdgeRemover'):

        super(SpuriousEdgeRemover, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID

        self.used_vertex_ids = set()
        self.found_homologies = list()


    def analyse(self):

        sortedVerts = sorted([x for x in self.graph.vertices], key=lambda x: len(self.graph.get_vertex(x).props['sequence']),
                             reverse=True)
        setRemoveVertexIDs = set()
        for x in sortedVerts:
            vertex = self.graph.get_vertex(x)
            targetEdges = []
            baseSeq = vertex.props['sequence']
            countArray = [0] * len(baseSeq)

            for edge in sorted(vertex.neighbors, key=lambda x: x.props['info'].identity, reverse=True):

                targetVertex = edge.target
                targetSeq = targetVertex.props['sequence']

                diamondResult = edge.props['info']
                targetVertexIDObj = self.getIDObj(edge, targetVertex)
                vertexIDObj = self.getIDObj(edge, vertex)

                if targetVertexIDObj != None:
                    for i in range(vertexIDObj.begin - 1, vertexIDObj.end):
                        countArray[i] += 1

            coverage = 0
            for x in countArray:
                if x > 0:
                    coverage += 1.0

            if vertex.name[1] == 'SE87_05730':
                pass
                # print(vertex)

            if coverage < len(baseSeq) / 2:
                # not enough coverage!
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

        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()


        return None