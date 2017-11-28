import copy

from utils import dictionaryHash, dictionaryEquals


class Edge:
    """

    Directed Edge class

    """
    def __init__(self, properties):

        self.props = {}

        if properties != None:
            self.props = properties

        self.source = None
        self.target = None

    def get_opposite_edge(self):
        """

        :return: shallow copy Edge
        """

        retEdge = Edge(self.props)
        retEdge.target = self.source
        retEdge.source = self.target

        return retEdge

    def getOpposite(self, elem):
        if elem == self.source:
            return self.target

        if elem == self.target:
            return self.source

        return None

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):

        if not isinstance(other, Edge):
            return False

        if self.source == other.source and self.target == other.target and self.props == other.props:
            return True

        return False

    def __hash__(self):
        """

        :return: hash for a directed edge (hash a->b != hash b->a
        """
        return self.source.__hash__() + 2*self.target.__hash__() + dictionaryHash(self.props)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.source.name) + " " + str(self.target.name) + " " + str(self.props)


class UndirectedEdge(Edge):
    """
    Undirected Edge class
    """
    def __eq__(self, other):
        if not isinstance(other, UndirectedEdge):
            return False

        if self.source == other.source and self.target == other.target and dictionaryEquals(self.props, other.props):
            return True

        if self.source == other.target and self.target == other.source and dictionaryEquals(self.props, other.props):
            return True

        return False

    def get_opposite_edge(self, copy_props=False):
        """

        :return: shallow copy Edge
        """

        props = self.props
        if copy_props:
            props = copy.deepcopy(self.props)

        retEdge = UndirectedEdge(props)
        retEdge.target = self.source
        retEdge.source = self.target

        return retEdge

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.source) + hash(self.target) + dictionaryHash(self.props)
