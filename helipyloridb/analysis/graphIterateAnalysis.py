from analysis.graphuser import GraphUser


class graphCleaner(GraphUser):

    def __init__(self, graph, genomeDB, stepID='graph_cleaner'):

        super(graphCleaner, self).__init__(graph, genomeDB, stepID)

        self.genomeDB = genomeDB
        self.stepID = stepID

    def analyse(self):
        self.graph.remove_empty_vertices()
        return None