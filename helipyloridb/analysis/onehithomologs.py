from analysis.graphuser import GraphUser
from analysis.homologyresults import HomologyResult

class OneHitHomologsConfig:

    def __init__(self, minIDScore=0.9, minLengthScore=0.8):

        self.minIDScore = minIDScore
        self.minLengthScore = minLengthScore



class oneHitHomologs(GraphUser):

    def __init__(self, graph, genomeDB, config, stepID='One2OneHomolog'):

        super(oneHitHomologs, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID

        self.config = config

        self.used_vertex_ids = set()
        self.found_homologies = list()

    def makeIdentityScore(self, alignment ):
        return alignment.identity

    def makeLengthScore(self, alignment ):

        qseq = self.genomeDB.get_sequence(alignment.query.genome, alignment.query.seqid)
        sseq = self.genomeDB.get_sequence(alignment.subject.genome, alignment.subject.seqid)

        lengthQuery =  (len(alignment.query) / len(qseq))
        lengthSubject =(len(alignment.subject) / len(sseq))

        return (lengthQuery+lengthSubject)/2.0

    def _analyse(self):

        retRes = HomologyResult()

        for vertexID in self.graph.vertices:

            if vertexID in self.used_vertex_ids:
                continue

            vertex = self.graph.vertices[vertexID]

            if vertex.name[1] == 'HP_0694':
                vertex.name = vertex.name

            if len(vertex.neighbors) == 1:

                targetVertex = vertex.neighbors[0].target

                # second condition is sanity check, should be always the case => unidirectional
                if len(targetVertex.neighbors) == 1 and targetVertex.neighbors[0].target == vertex:
                    diamondResult = targetVertex.neighbors[0].props['info']
                    identityScore = self.makeIdentityScore(diamondResult) > self.config.minIDScore
                    lengthScore = self.makeLengthScore(diamondResult) > self.config.minLengthScore

                    if identityScore and lengthScore:
                        self.used_vertex_ids.add(vertex.name)
                        self.used_vertex_ids.add(targetVertex.name)

                        #print("Step2", vertex.name, targetVertex.name, diamondResult)

                        retRes.homology_relations.append( (vertex.name, targetVertex.name,
                                                        {'step': self.stepID, 'file': "n/a", 'edge': (vertex.name, targetVertex.name)}
                                                        )
                                                      )

        self.graph.remove_vertices_by_id(self.used_vertex_ids)
        self.graph.remove_empty_vertices()

        return retRes