import io
from collections import defaultdict

from utils import mergeDicts


class MatchingRegion:

    def __init__(self, range1, range2):

        self.range1 = range1
        self.range2 = range2



class MultiCombination:

    def __init__(self, props=None):

        self.props = props
        self.elems = []

    def addMatch(self, elemInterval1, elemInterval2):

        match = MatchingRegion(elemInterval1, elemInterval2)
        self.elems.append(match)


class HomologyDatabase:

    def __init__(self):

        self.homologies = dict()
        self.homologyProperties = defaultdict(lambda: dict())

        self.combinations = defaultdict(set)
        self.combinationProperties = defaultdict(lambda: dict())

        self.multiCombinations = []

    def findHomologyForID(self, searchID):

        for homolID in self.homologies:

            elems = self.homologies[homolID]

            for id in elems:
                if id == searchID:
                    return homolID

        return None

    def addMultiCombination(self, multiCombi):

        if not isinstance(multiCombi, MultiCombination):
            raise ValueError("multiCombi must be MultiCombination")

        self.multiCombinations.append(multiCombi)

    def addHomologyRelation(self, id1, id2, properties=None):

        if id1 < id2:
            idtuple = (id1, id2)
        else:
            idtuple = (id2, id1)

        for x in self.homologies:

            allElems = self.homologies[x]

            if id1 in allElems or id2 in allElems:
                allElems.add(id1)
                allElems.add(id2)

                self.homologies[x] = allElems
                if properties != None:
                    self.homologyProperties[x][ idtuple ] = mergeDicts(self.homologyProperties[x].get(idtuple, None), properties)
                return

        newRel = set()
        newRel.add(id1)
        newRel.add(id2)

        homID = "HOMID" + str(len(self.homologies))
        self.homologies[homID] = newRel
        if properties != None:
            self.homologyProperties[homID][ idtuple ] = properties

    def finalize(self):

        changed = True
        while changed:

            allIDs = [x for x in self.homologies]
            changed = False

            for i in range(0, len(allIDs)):
                if changed:
                    break

                for j in range(i+1, len(allIDs)):

                    xi = allIDs[i]
                    xj = allIDs[j]

                    setGenesI = self.homologies[xi]
                    setGenesJ = self.homologies[xj]

                    if len(setGenesI.intersection(setGenesJ)) > 0:
                        newSet = setGenesI.union(setGenesJ)
                        self.homologies[xi] = newSet
                        del self.homologies[xj]
                        changed=True

                        break




    def addCombination(self, id1, listIDs, properties):

        for x in listIDs:
            self.combinations[id1].add(x)

        self.combinationProperties[id1] = mergeDicts(self.combinationProperties[id1], properties)

    def printCombinations(self):

        for x in self.combinations:

            for partner in self.combinations[x]:
                print("\t".join([str(xn) for xn in x]) + "\t" + "\t".join([str(xp) for xp in partner]))

    def __str__(self):


        outStr = io.StringIO()

        for homid in self.homologies:

            allRels = self.homologies[homid]
            allProps = self.homologyProperties[homid]

            propsAssigned = None

            for rel in allRels:

                if allProps != None and len(allProps) > 0:
                    propKeys = [x for x in allProps]

                    for propKey in propKeys:
                        if rel in propKey:
                            propsAssigned = allProps[propKey]
                            break

                attribsPrint = [str(x) for x in rel]
                if propsAssigned != None:
                    attribsPrint.append(str(propsAssigned))


                outStr.write(homid + "\t" + "\t".join(attribsPrint) + "\n")


        outString = outStr.getvalue()
        outStr.close()

        return outString


    @classmethod
    def loadFromFile(cls, filename):

        homdb = HomologyDatabase()
        homdb.homologies = defaultdict(set)

        with open(filename, 'r') as infile:

            for line in infile:

                if len(line) == 0 or len(line.strip()) == 0:
                    continue

                aline = line.strip().split('\t')

                homdb.homologies[aline[0]].add( (aline[1], aline[2]) )

        return homdb