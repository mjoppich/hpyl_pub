import codecs
from collections import OrderedDict, Counter


import csv

from io import StringIO

import re


class Synonym:

    @classmethod
    def parseFromLine(cls, line):
        aLine = line.split(':', 1)

        retSyn = Synonym(aLine[0].strip())

        aSyns = aLine[1].split('|')
        for x in [x.strip() for x in aSyns]:
            retSyn.syns.append(x)

        return retSyn

    def __init__(self, id):

        id = id.strip()
        id = id.replace(':', '_')
        self.id = id
        self.currentIdx = 0
        self.syns = []

    def get(self, synIdx):

        if synIdx >= 0 and synIdx < len(self.syns):
            return self.syns[synIdx]

        return None



    def removeCommonSynonymes(self, commonSyns):

        for x in commonSyns:
            if x in self:
                self.removeSyn(x)

    def getSynonymes(self):

        return self.syns

    def addTextSyns(self, synText):

        if synText == None or synText == 'None':
            return

        synText = synText.strip()

        if len(synText) == 0:
            return

        vAllWords = self.getAllSplittedSyns(synText, ', ')

        for x in vAllWords:
            self.addSyn(x)


    def addSyn(self, newSyn):

        if newSyn == None:
            return

        newSyn = newSyn.strip()

        if len(newSyn) == 0:
            return
        if len(newSyn) < 3:
            return

        if newSyn[0] == newSyn[len(newSyn)-1] and newSyn[0] == '"':
            newSyn = newSyn[1:len(newSyn)-1]

        if len(newSyn) < 3:
            return
        if len(newSyn) == 0:
            return

        if newSyn.startswith('symbol withdrawn, '):
            newSyn = newSyn[len('symbol withdrawn, '):]

        if newSyn.startswith("see "):
            newSyn = newSyn[len("see "):]

        if newSyn.endswith("~withdrawn"):
            newSyn = newSyn[:-len("~withdrawn")]

        if '|' in newSyn:
            newSyn = newSyn.replace('|', '-')

        if newSyn.upper() in ['SYMBOL WITHDRAWN', 'WITHDRAWN', 'PSEUDOGENE', 'RNA', 'ENTRY WITHDRAWN']:
            return

        newSyn = re.sub('\s\s+',' ',newSyn)

        if len(re.sub('[0-9\W]*', '', newSyn)) == 0:
            return

        if newSyn == '1,3':
            print(newSyn)

        if not newSyn in self.syns:
            self.syns.append(newSyn)

    def __contains__(self, item):
        return item in self.syns

    def __iter__(self):

        self.currentIdx = 0
        return self

    def __next__(self):

        if self.currentIdx >= len(self.syns):
            raise StopIteration

        self.currentIdx+=1
        return self.syns[self.currentIdx-1]

    def __str__(self):
        return self.id.replace(':', '_') + ":" + "|".join(self.syns)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.syns)

    def removeSynUpper(self, dictOfSyns, removeExtraChars=False):

        for excludeName in dictOfSyns:

            listToExclude = dictOfSyns[excludeName]

            syns2remove = set()
            for syn in self.syns:

                usyn = syn.upper()

                if syn in listToExclude or usyn in listToExclude:
                    syns2remove.add(syn)

            for syn in syns2remove:
                self.removeSyn(syn)

    def addAlphaBetaVariants(self):
        self.__searchReplaceVariant('A', 'alpha', [u"\u03B1", '-'+u"\u03B1"])
        self.__searchReplaceVariant('B', 'beta', [u"\u03B2", '-' + u"\u03B2"])


    def __searchReplaceVariant(self, endChar, endWord, replaceWith):

        endsWithChar = self.id.endswith(endChar)
        endsWithWord = False

        for x in self.syns:
            if x.endswith(endWord) or x.upper().endswith(endWord.upper()):
                endsWithWord = True
                break

        if endsWithChar and endsWithWord:

            for x in replaceWith:
                aid = self.id.rsplit(endChar, 1)
                aid.append(x)

                newsyn = "".join(aid)

                print(self.id + " ADD SYN " + newsyn)
                self.syns.append(newsyn)


    def removeNumbers(self):

        i = 0
        while i < len(self.syns):

            try:
                float(self.syns[i])
                self.removeSyn(i)

            finally:
                i+= 1


    def removeSyn(self, syn):

        if type(syn) == str:

            while syn in self.syns:
                idx = self.syns.index(syn)
                self.removeSyn(idx)

        else:

            isyn = int(syn)

            if isyn >= 0 and isyn < len(self.syns):
                del self.syns[isyn]

    def splitQuotedDelimited(self, search, delimiter=', ', quotechars=['\"', '\'']):

            allWords = []
            i = 0
            while i <len(search):

                x = search[i]

                if x in quotechars:
                    y = search.find(x, i+1)

                    if y != -1:
                        allWords.append((i,y))
                        i = y

                i += 1


            foundWords = []

            i = -len(delimiter)
            lasti = 0
            lastAdded = 0
            while i < len(search):

                lasti = i+len(delimiter)
                i = search.find(delimiter, lasti)

                if i == -1:
                    i = len(search)
                    break

                if len(allWords) > 0:
                    ignoreDel = False
                    for word in allWords:
                        if word[0] <= i and i <= word[1]:
                            ignoreDel = True
                            break

                    if ignoreDel:
                        continue

                foundWords.append(search[lastAdded:i])
                lastAdded = i+len(delimiter)

            if lastAdded != len(search):
                foundWords.append(search[lastAdded:len(search)])

            return foundWords


    def findNestedMatches(self, search, delimiter=',', quotechars={'(': ')'}):

        nestLevel = 0
        nestStart = 0
        stack = []

        foundNesting = []

        i = 0
        while i < len(search):

            char = search[i]

            if nestLevel == 0:
                nestStart = i

            if char in quotechars:
                nestLevel += 1
                stack.append( quotechars[char] )

            elif len(stack) > 0 and char == stack[-1]:
                stack = stack[0:len(stack)-1]
                nestLevel -= 1

                if nestLevel == 0:

                    nestWord = search[nestStart+1:i]
                    foundNesting.append(nestWord)

            i += 1

        return foundNesting



    def getAllSplittedSyns(self, search, delimiter=', ', quotechars=['\"', '\'']):

        vAllWords = self.splitQuotedDelimited(search)

        setAllWords = set()

        for word in vAllWords:
            setAllWords.add(word)

            foundBrackets = []

            foundBracketWords = self.findNestedMatches(word)
            for x in foundBracketWords:

                foundBrackets.append( "("+x+")" )

                for y in self.splitQuotedDelimited(x):

                    if y[0]==y[len(y)-1] and len(y) > 0:
                        y = y[1:len(y)-1]

                    setAllWords.add(y)

            testWord = word
            for bracketWord in foundBrackets:
                testWord = testWord.replace(bracketWord, '')

            setAllWords.add(testWord)

        return list(setAllWords)

class Synfile:

    def __init__(self, sFileLocation):

        self.mSyns = {}
        self.line2syn = {}

        def addSyn(sLine, iLine):

            oSyn = Synonym.parseFromLine(sLine)

            self.mSyns[ oSyn.id ] = oSyn
            self.line2syn[iLine] = oSyn.id

        with codecs.open(sFileLocation, 'r', 'latin1') as infile:
            idx = 0
            for line in infile:
                addSyn(line, idx)
                idx += 1

        self.synIDs = None
        self.synIDidx = None

    def __iter__(self):

        self.synIDs = [x for x in self.mSyns]
        self.synIDidx = 0

        return self

    def __next__(self):

        curIdx = self.synIDidx
        self.synIDidx += 1

        if curIdx < len(self.synIDs):
            return self.mSyns[self.synIDs[curIdx]]

        raise StopIteration()


    def __len__(self):
        return len(self.mSyns)

    def get(self, iSynID):

        return self.mSyns.get(self.line2syn.get(iSynID, None), None)

    def histogramSynonymes(self):

        oSynCounter = Counter()

        for iSynID in self.mSyns:

            oSyn = self.mSyns[iSynID]

            vSyns = set(oSyn.getSynonymes())

            for sSyn in vSyns:

                oSynCounter[sSyn] += 1

        return oSynCounter

    def histogramSynonymesOrdered(self, limit = None):

        synCounter = self.histogramSynonymes()

        if limit != None:
            synOrdered = OrderedDict( synCounter.most_common(limit) )
        else:
            synOrdered = OrderedDict(sorted(synCounter.items(), key=lambda x: x[1]))

        return synOrdered