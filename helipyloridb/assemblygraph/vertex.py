from assemblygraph.edge import Edge, UndirectedEdge


class Vertex:
    def __init__(self, vertex, properties):
        self.name = vertex
        self.props = properties
        self.neighbors = []

    def __eq__(self, other):

        if hash(self) != hash(other):
            return False

        if not isinstance(other, Vertex):
            return False

        return other.name == self.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.name)# + dictionaryHash(self.props)# + hash(tuple(self.neighbors))

    def connectsTo(self, other):

        if isinstance(other, Vertex):
            othersNames = [other.name]
        elif isinstance(other, (list, tuple, set)):
            othersNames = other

        for edge in self.neighbors:
            if edge.target.name in othersNames:
                return True
            elif edge.source.name in othersNames:
                return True

        return False

    def add_neighbor(self, neighbor, props=None, undirected=False):
        if isinstance(neighbor, Vertex):

            if undirected:
                edge = UndirectedEdge(props)
            else:
                edge = Edge(props)

            edge.source = self
            edge.target = neighbor
            self.neighbors.append(edge)
            return edge

        elif isinstance(neighbor, Edge) and props == None:
            neighbor.source = self
            self.neighbors.append(neighbor)
            return neighbor

        raise ValueError("neighbor must be either Vertex or Edge")

    def __str__(self):
        return str(self.name) + " " + str(self.props) + " " + str(self.neighbors)

    def __repr__(self):
        return self.__str__()