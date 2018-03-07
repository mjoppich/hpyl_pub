import logging
import os
from collections import Counter

import editdistance
import re

from xrefs.TextmineDocument import TextMineDocument
from xrefs.tmSentences import SentenceID, Sentence

logger = logging.getLogger('convertJatsToText')

import nltk.data


class Sentenicer:

    def __init__(self, tokenizer_loc = 'tokenizers/punkt/english.pickle'):

        self.tokenizer = nltk.data.load(tokenizer_loc)
        self.ignoredSections = Counter()
        pass


    def _splitToSentences(self, content):


        #aPseudoSentences = content.split(".")
        #vSents = [x.strip() for x in aPseudoSentences if len(x) > 7]

        vSents = self.tokenizer.tokenize(content)

        return vSents

    def _makeSentences(self, articleName, module, sents):

        iSent = 0
        allSentences = []

        for x in sents:

            x = x.strip()
            x = x.strip(',.;')

            if len(x) > 0:

                sentID = SentenceID(articleName, module, iSent)
                thisSentence = Sentence(sentID, x)

                allSentences.append(thisSentence)

                iSent += 1

        return allSentences

    def editDist(self, x, y):
        return editdistance.eval(x.upper(),y.upper())


    def run(self, doc, ignoreSections = list(), includeReferences=True):

        if not isinstance(doc, TextMineDocument):
            raise TypeError("doc must be TextMineDocument" + str(doc))

        ignoreSections = [x.upper() for x in ignoreSections]

        allCreatedSentences = []

        titles = " ".join(doc.titles())
        abstracts = " ".join(doc.abstracts())

        titleSents = self._splitToSentences(titles)
        abstractSents = self._splitToSentences(abstracts)

        if (len(titleSents) == 0):
            logger.info("Titles empty for article: " + str(doc.id()))

        if (len(abstractSents) == 0):
            logger.info("Abstracts empty for article: " + str(doc.id()))

        if len(titleSents) > 0:
            allCreatedSentences += self._makeSentences(doc.name(), 1, titleSents)

        if len(abstractSents) > 0:
            allCreatedSentences += self._makeSentences(doc.name(), 2, abstractSents)


        ### PREPARE CONTENT
        docSections = doc.sections()

        contentList = []
        if len(docSections) > 0:
            ## check which sections to take
            for x in docSections:

                sectionTitle = re.sub(r'[\W_]+', '', x)

                addSection = True
                for sect in ignoreSections:

                    if self.editDist(sectionTitle.upper(), sect) <= 3:
                        addSection = False
                        self.ignoredSections[sect] += 1
                        break

                if addSection:
                    contentList.append(docSections[x])

        else: #if len(contentList) == 0:
            contentList = doc.texts()

        content = " ".join(contentList)
        contentSents = self._splitToSentences(content)

        if (len(contentSents) == 0):
            logger.info("Contents empty for article: " + str(doc.id()))

        if len(contentSents) > 0:
            allCreatedSentences += self._makeSentences(doc.name(), 3, contentSents)

        if includeReferences:
            references = " ".join(doc.reflists())
            refSents = self._splitToSentences(references)
            allCreatedSentences += self._makeSentences(doc.name(), 4, refSents)

        return allCreatedSentences

    def printIgnoration(self):

        for (section, cnt) in self.ignoredSections.most_common():
            print(str(section) + " " + str(cnt))