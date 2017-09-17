class Edge:

    def __init__(self, properties):

        self.props = {}

        if properties != None:
            self.props = properties

        self.source = None
        self.target = None
