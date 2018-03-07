from collections import defaultdict


class SentenceID:

    def __init__(self, docID, paraNum, sentNum):

        self.docid = docID
        self.paragraph = int(paraNum)
        self.sentnum = int(sentNum)

    def __str__(self):

        return str(self.docid) + "." + str(self.paragraph) + "." + str(self.sentnum)

    def __repr__(self):

        return str(self)

    @classmethod
    def fromString(cls, input):

        if type(input) == str:
            aids = input.split('.')
        else:
            aids = input

        return SentenceID(aids[0], int(aids[1]), int(aids[2]))

class Sentence:

    def __init__(self, sentID, sent):

        if not isinstance(sentID, SentenceID):
            raise ValueError("sentID must be SentenceID: " + str(sentID))

        self.id = sentID
        self.text = sent

    def __str__(self):

        return str(self.id) + "\t" + str(self.text)

    def __repr__(self):
        return self.__str__()

    @classmethod
    def fromLine(cls, line):
        aline = line.split('\t')

        sentID = SentenceID(aline[0])
        sent = aline[1]

        return Sentence(sentID, sent)


class SentencesFile:

    def __init__(self, sFileLocation):

        self.mSents = defaultdict(list)

        def addSent(sLine, iLine):

            oSent = Sentence(sLine)
            self.mSents.addVectorized(oSent.id.docid, oSent)

        fu.readFile(sFileName=sFileLocation, sFunc=addSent)

    def get(self, docid):

        if docid in self.mSents:

            return self.mSents[docid]

        return None

    def __len__(self):

        length = 0
        for x in self.mSents:
            length += len(self.mSents[x])

        return len


    def getDocumentCount(self):

        return len(self.mSents)
