from analysis.graphuser import GraphUser
from analysis.homologyresults import HomologyResult


class oneMultipleConfig:

    def __init__(self):
        self.minIdentity = 0.8
        self.minQueryLength = 0.85
        self.minSubjectLength = 0.85
        self.edgeSortExpression = lambda x: x.props['info'].identity
        self.allowPartialLength = False
        self.allowMultiple = False
        self.betterEdgeCheck = False

class oneMultipleHomologs(GraphUser):

    def __init__(self, graph, genomeDB, config, stepID='one2multiHitHomologs'):

        super(oneMultipleHomologs, self).__init__(graph, genomeDB, stepID)
        self.config = config
        self.stepID = stepID

        self.used_vertex_ids = set()
        self.found_homologies = list()

    def makeIdentityScore(self, alignment ):
        return alignment.identity

    def makeLengthScore(self, alignment ):

        qseq = self.genomeDB.get_sequence(alignment.query.genome, alignment.query.seqid)
        sseq = self.genomeDB.get_sequence(alignment.subject.genome, alignment.subject.seqid)

        lengthQuery =  (len(alignment.query) / len(qseq))
        lengthSubject =(len(alignment.subject) / len(sseq))

        if lengthQuery < 0.8:
            return 0
        if lengthSubject < 0.8:
            return 0

        return (lengthQuery+lengthSubject)/2.0

    def _analyse(self):

        returnResult = HomologyResult()

        sortedVerts = sorted([x for x in self.graph.vertices], key=lambda x: len(self.graph.get_vertex(x).props['sequence']),
                             reverse=True)
        setRemoveVertexIDs = set()

        for x in sortedVerts:
            vertex = self.graph.get_vertex(x)

            if vertex.name in setRemoveVertexIDs:
                continue # already assigned

            if vertex.name[1] == 'HPP12_1154':
                vertex.name = vertex.name

            for edge in sorted(vertex.neighbors, key=self.config.edgeSortExpression, reverse=True):

                targetVertex = edge.target

                if targetVertex.name in setRemoveVertexIDs:
                    continue # already assigned

                vertexIDObj = self.getIDObj(edge, vertex)
                targetVertexIDObj = self.getIDObj(edge, targetVertex)

                queryLength = len(vertexIDObj) / len(vertex.props['sequence'])
                subjectLength = len(targetVertexIDObj) / len(targetVertex.props['sequence'])

                acceptEdge = False

                considerEdge = False
                if self.config.allowPartialLength:
                    considerEdge = (queryLength > self.config.minQueryLength or subjectLength > self.config.minSubjectLength) and edge.props[
                                                                                                              'info'].identity > self.config.minIdentity
                    considerEdge = considerEdge and min([queryLength, subjectLength]) > 0.5
                else:
                    considerEdge = queryLength > self.config.minQueryLength and subjectLength > self.config.minSubjectLength and edge.props[
                                                                                                             'info'].identity > self.config.minIdentity

                if considerEdge:

                    acceptEdge = True

                    for targetEdge in targetVertex.neighbors:
                        if targetEdge.props['info'].identity > edge.props['info'].identity:
                            otherVertexObj = self.getIDObj(targetEdge, targetEdge.target)
                            otherVertexLength = len(otherVertexObj) / len(targetEdge.target.props['sequence'])

                            if subjectLength > 0.9 and otherVertexLength > 0.9:
                                continue

                            if otherVertexLength > subjectLength:
                                acceptEdge = False
                                break

                if acceptEdge:

                    if vertex.name in setRemoveVertexIDs or targetVertex.name in setRemoveVertexIDs:
                        self.log_warn("Duplicate vertices: " + str(vertex.name))
                        self.log_warn("Duplicate vertices: " + str(targetVertex.name))

                    setRemoveVertexIDs.add(vertex.name)
                    setRemoveVertexIDs.add(targetVertex.name)
                    # print("acceptOneOfMultiple", minIdentity, minQueryLength, minSubjectLength, allowPartialLength, vertex.name, targetVertex.name, edge.props['info'])

                    if vertex.name[1] == 'HP_0694':
                        vertex.name = vertex.name

                    returnResult.homology_relations.append( (vertex.name, targetVertex.name,
                                                {'step': self.stepID, 'file': "N/A",
                                                 'edge': (vertex.name, targetVertex.name)}
                                                ) )

                    if not self.config.allowMultiple:
                        break

        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()

        return returnResult