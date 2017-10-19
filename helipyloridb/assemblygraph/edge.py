class Edge:

    def __init__(self, properties):

        self.props = {}

        if properties != None:
            self.props = properties

        self.source = None
        self.target = None

    def __str__(self):

        return str(self.source) + " " + str(self.target) + " " + str(self.props)
