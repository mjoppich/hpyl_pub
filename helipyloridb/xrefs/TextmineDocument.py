from abc import ABC, abstractmethod

class TextMineDocument(ABC):

    @abstractmethod
    def abstracts(self):
        """

        :return: list of abstracts for given document
        """
        pass

    @abstractmethod
    def titles(self):
        """

        :return: list of titles for given document
        """
        pass

    @abstractmethod
    def id(self):
        pass

    @abstractmethod
    def name(self):
        pass

    @abstractmethod
    def sections(self):
        """

        :return: list of sections for given document
        """
        pass

    @abstractmethod
    def reflists(self):
        """

        :return: list of ref elemnts/literature for given document
        """
        pass

    @abstractmethod
    def texts(self):
        """

        :return: list of texts for given document
        """
        pass
