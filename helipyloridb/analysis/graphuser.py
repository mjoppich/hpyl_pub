

from abc import ABC, abstractmethod

import logging


class GraphUser(ABC):

    def __init__(self, graph, genomeDB, stepID = "GraphUser"):

        self.graph = graph
        self.stepID = stepID
        self.genomeDB = genomeDB

        self.logger = None

        self._init_logger()

    def _init_logger(self):

        FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(format=FORMAT)

        self.logger = logging.getLogger( self.stepID )


    def log_warn(self, message):
        self.logger.warn(message)
        
    def log_info(self, message):
        self.logger.info(message)


    @abstractmethod
    def analyse(self):
        """

        :return: returns objects for homology database
        """
        pass

    def getIDObj(self, edge, vertex):

        diamondResult = edge.props['info']

        if vertex.name == (diamondResult.query.genome, diamondResult.query.seqid):
            return diamondResult.query

        if vertex.name == (diamondResult.subject.genome, diamondResult.subject.seqid):
            return diamondResult.subject

        return None

