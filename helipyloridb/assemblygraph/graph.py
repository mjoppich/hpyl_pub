from assemblygraph.vertex import Vertex


class Graph:
    def __init__(self):
        self.vertices = {}

    def get_vertices(self):
        return [x for x in self.vertices]

    def add_vertex(self, vertex, replace=False):
        if isinstance(vertex, Vertex):

            if replace == False and vertex.name in self.vertices:
                raise ValueError("Node already exists: " + vertex.name)

            self.vertices[vertex.name] = vertex

            return vertex

        raise ValueError("vertex must be a Vertex")

    def add_vertex_if_not_exists(self, vertex):

        if not self.contains_vertex(vertex):
            self.add_vertex(vertex)
            return vertex
        else:
            return self.get_vertex(vertex)

    def add_vertices(self, vertices, replace=False):
        for vertex in vertices:
            self.add_vertex(vertex, replace)

    def contains_vertex(self, vertex_id):

        if self.vertices.get(vertex_id, None) != None:
            return True

        if isinstance(vertex_id, Vertex):
            return vertex_id.name in self.vertices
        #else:
        #    return vertex_id in self.vertices

    def get_vertex(self, vertex_id, default=None):

        if isinstance(vertex_id, Vertex):
            return self.vertices.get(vertex_id.name, default)
        else:
            return self.vertices.get(vertex_id, default)

    def add_edge(self, vertex_from, vertex_to, edge_props, undirected=False):

        if not isinstance(vertex_from, Vertex) or not isinstance(vertex_to, Vertex):
            raise ValueError("Both nodes must be vertices")

        if not self.contains_vertex(vertex_from):
            self.add_vertex(vertex_from)
        if not self.contains_vertex(vertex_to):
            self.add_vertex(vertex_to)

        retEdge = vertex_from.add_neighbor(vertex_to, edge_props, undirected=undirected)

        if undirected:
            retOppEdge = retEdge.get_opposite_edge()
            vertex_to.add_neighbor(retOppEdge)

        return retEdge

    def get_edges(self, decision_function=lambda edge: True):

        allEdges = set()

        for vertex in self.vertices:
            vertexNode = self.vertices[vertex]

            for edge in vertexNode.neighbors:

                if decision_function(edge):
                    allEdges.add(edge)

        return allEdges

    def remove_vertices_by_id(self, vertixIDs):

        for vid in self.vertices:

            if vid in vertixIDs:
                continue

            vertNode = self.vertices[vid]

            if not vertNode.connectsTo( vertixIDs ):
                continue

            remEdgesIdx = set()
            for edge in vertNode.neighbors:
                if edge.target.name in vertixIDs:
                    remEdgesIdx.add(edge)
                elif edge.source.name in vertixIDs:
                    remEdgesIdx.add(edge)

            for x in remEdgesIdx:
                idx = vertNode.neighbors.index(x)
                del vertNode.neighbors[idx]

        for vid in vertixIDs:
            if vid in self.vertices:
                del self.vertices[vid]


    def remove_vertex(self, vertex):

        if self.contains_vertex(vertex):
            vertexObj = self.get_vertex(vertex)

            for vid in self.vertices:

                if vid == vertex:
                    continue

                vert = self.vertices[vid]

                if not vert.connectsTo(vertexObj):
                    continue

                remEdgesIdx = set()
                for edge in vert.neighbors:
                    if edge.target.name == vertexObj.name:
                        remEdgesIdx.add(edge)
                    elif edge.source.name == vertexObj.name:
                        remEdgesIdx.add(edge)

                for x in remEdgesIdx:
                    idx = vert.neighbors.index(x)
                    del vert.neighbors[idx]

            del self.vertices[vertexObj.name]

            return

        raise ValueError("Vertex does not exist", vertex)

    def _propsMatch(self, lprops, rprops):

        for x in lprops:
            if not x in rprops:
                return False

            if lprops[x] != rprops[x]:
                return False

        return True

    def vertex_count(self):

        return len(self.vertices)

    def edge_count(self):

        cnt = 0
        for vertexID in self.vertices:

            cnt += len(self.vertices[vertexID].neighbors)

        return cnt

    def remove_empty_vertices(self):

        removeVertexIDs = set()
        for vertexID in self.vertices:

            vertex = self.vertices[vertexID]

            if len(vertex.neighbors) == 0:
                removeVertexIDs.add(vertexID)

        self.remove_vertices_by_id(removeVertexIDs)


    def remove_edge(self, vertex_from, vertex_to, edge_props, undirected=False):

        if not isinstance(vertex_from, Vertex) or not isinstance(vertex_to, Vertex):
            raise ValueError("Both nodes must be vertices")

        if not self.contains_vertex(vertex_from) or not self.contains_vertex(vertex_to):
            raise ValueError("Both nodes must already exist")

        setRemoveEdges = set()
        for edge in vertex_from.neighbors:
            if edge.target == vertex_to:

                if self._propsMatch(edge_props, edge.props):
                    setRemoveEdges.add(edge)

        for edge in setRemoveEdges:
            vertex_from.neighbors.remove(edge)

        if undirected==True:
            self.remove_edge(vertex_to, vertex_from, edge_props, False)

    def add_edges(self, edges):
        for edge in edges:
            self.add_edge(edge.source, edge.target, edge_props=edge.props)
