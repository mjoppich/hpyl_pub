from analysis.graphuser import GraphUser, IDUser
from analysis.homologyresults import HomologyResult


class oneMultipleConfig(IDUser):

    def __init__(self):
        self.minIdentity = 0.8
        self.minQueryLength = 0.85
        self.minSubjectLength = 0.85
        self.edgeSortExpression = lambda x: x.props['info'].identity
        self.allowPartialLength = False
        self.allowMultiple = False
        self.betterEdgeCheck = False

        def defaultEdgeFunc(confObj, edge, source, target):
            queryLength = confObj.get_seq_fraction(edge, source)
            subjectLength = confObj.get_seq_fraction(edge, target)

            retVal = queryLength > confObj.minQueryLength#
            retVal = retVal and subjectLength > confObj.minSubjectLength
            retVal = retVal and edge.props['info'].identity > confObj.minIdentity

            return retVal

        self.considerEdgeFunc = defaultEdgeFunc


    def considerEdge(self, edge, source, target):

        return self.considerEdgeFunc(self, edge, source, target)


class oneMultipleHomologs(GraphUser):

    def __init__(self, graph, genomeDB, config, stepID='OneOfManyHomolog'):

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

            if not self.config.allowMultiple and vertex.name in setRemoveVertexIDs:
                continue # already assigned

            if vertex.name[1] in ['HP_1526', 'jhp_0476']:
                vertex.name = vertex.name

            for edge in sorted(vertex.neighbors, key=self.config.edgeSortExpression, reverse=True):

                targetVertex = edge.target

                if not self.config.allowMultiple and targetVertex.name in setRemoveVertexIDs:
                    continue # already assigned

                acceptEdge = False

                considerEdge = self.config.considerEdge(edge, vertex, targetVertex)
                subjectLength = self.get_seq_fraction(edge, targetVertex)

                if considerEdge:

                    acceptEdge = True

                    for targetEdge in targetVertex.neighbors:
                        if targetEdge.props['info'].identity > edge.props['info'].identity:

                            otherVertexLength = self.get_seq_fraction(targetEdge, targetEdge.target)

                            if subjectLength > self.config.minQueryLength and otherVertexLength > self.config.minSubjectLength:
                                continue

                            if otherVertexLength > subjectLength:
                                acceptEdge = False
                                break

                if acceptEdge:

                    if not self.config.allowMultiple and (vertex.name in setRemoveVertexIDs or targetVertex.name in setRemoveVertexIDs):
                        if vertex.name in setRemoveVertexIDs:
                            self.log_warn("Duplicate vertices: " + str(vertex.name))
                            self.log_warn("Opposite: " + str(targetVertex.name))
                        elif targetVertex.name in setRemoveVertexIDs:
                            self.log_warn("Duplicate vertices: " + str(targetVertex.name))
                            self.log_warn("Opposite: " + str(vertex.name))
                        else:
                            self.log_warn("Duplicate vertex: " + str(targetVertex.name))
                            self.log_warn("Duplicate vertex: " + str(vertex.name))

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
                else:
                    pass

        self.graph.remove_vertices_by_id(setRemoveVertexIDs)
        self.graph.remove_empty_vertices()

        return returnResult