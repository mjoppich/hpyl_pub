
class Edge:

    def __init__(self, properties):

        self.props = {}

        if properties != None:
            self.props = properties

        self.source = None
        self.target = None


    def getOpposite(self, elem):
        if elem == self.source:
            return self.target

        if elem == self.target:
            return self.source

        return None


    def __str__(self):
        return str(self.source) + " " + str(self.target) + " " + str(self.props)


class UndirectedEdge(Edge):
    def __hash__(self):
        return self.source.__hash__() + self.target.__hash__()