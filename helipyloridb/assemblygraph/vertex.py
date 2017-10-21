from assemblygraph.edge import Edge, UndirectedEdge


class Vertex:
    def __init__(self, vertex, properties):
        self.name = vertex
        self.props = properties
        self.neighbors = []

    def add_neighbor(self, neighbor, props=None, undirected=False):
        if isinstance(neighbor, Vertex):

            if undirected:
                edge = UndirectedEdge(props)
            else:
                edge = Edge(props)

            edge.source = self
            edge.target = neighbor
            self.neighbors.append(edge)
            return

        elif isinstance(neighbor, Edge) and props == None:
            neighbor.source = self
            self.neighbors.append(neighbor)
            return

        raise ValueError("neighbor must be either Vertex or Edge")

    def __str__(self):
        return str(self.name) + " " + str(self.props) + " " + str(self.neighbors)

    def __repr__(self):
        return self.__str__()