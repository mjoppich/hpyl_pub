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

    def to_db_str(self, combID):

        outstr=""
        for match in self.elems:
            outstr += "{id}\t{gen}\t{nam}\n".format(id=combID, gen=match.range1.genome, nam=match.range1.seqid)
            outstr += "{id}\t{gen}\t{nam}\n".format(id=combID, gen=match.range2.genome, nam=match.range2.seqid)

        return outstr


class HomologyDatabase:

    def __init__(self):

        self.homologies = dict()
        self.homologyProperties = defaultdict(lambda: dict())

        self.combinations = defaultdict(set)
        self.combinationProperties = defaultdict(lambda: dict())

        self.multiCombinations = []

    def get_all_ids(self):

        return [x for x in self.homologies] + [x for x in self.combinations]# + [x for x in self.multiCombinations]

    def get_all_organisms(self):

        orgs = set()

        for x in self.homologies:
            for y in self.homologies[x]:
                orgs.add(y[0])

        return orgs


    def get_cluster(self, homID):

        clusterProt = self.homologies.get(homID, None)

        if clusterProt == None:
            clusterProt = self.combinations.get(homID, None)

        if clusterProt == None:
            return None

        clusterDict =  defaultdict(list)

        for elem in clusterProt:
            clusterDict[elem[0]].append(elem[1])

        return clusterDict

    def get_homology_cluster(self, homID):

        clusterProt = self.homologies.get(homID)

        if clusterProt == None:
            return None

        clusterDict =  {}

        for elem in clusterProt:
            clusterDict[elem[0]] = elem[1]

        return clusterDict

    def get_organism_elements(self, orgID):

        foundElems = set()

        for homolID in self.homologies:
            elems = self.homologies[homolID]

            for id in elems:
                if id[0] == orgID:
                    foundElems.add(id[1])

        for homolID in self.combinations:
            elems = self.combinations[homolID]
            for id in elems:
                if id[0] == orgID:
                    foundElems.add(id[1])

        return foundElems

    def get_hom_ids(self):
        return [x for x in self.homologies]

    def get_comb_ids(self):
        return [x for x in self.combinations]

    def get_mulcombs(self):
        return [x for x in self.multiCombinations]

    def findHomologyForGeneID(self, searchID):

        retHomIDs = []

        for homolID in self.homologies:

            elems = self.homologies[homolID]

            for id in elems:
                if id[1] == searchID:
                    retHomIDs.append( homolID )

        return list(set(retHomIDs))

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

                if len(allElems) > 100:
                    print(id1, id2)

                self.homologies[x] = allElems
                if properties != None:
                    self.homologyProperties[x][ idtuple ] = mergeDicts(self.homologyProperties[x].get(idtuple, None), properties)

                return

        newRel = set()
        newRel.add(id1)
        newRel.add(id2)


        if len (self.homologies) > 0:
            allIDs = [int(x.replace('HOMID', '')) for x in self.homologies]
            maxID = max(allIDs) + 1
        else:
            maxID = 1

        homID = "HOMID" + str(maxID)
        self.homologies[homID] = newRel
        if properties != None:
            self.homologyProperties[homID][ idtuple ] = properties

    def finalize(self):

        #fetch all combid's
        #for combElem in self.combinations:
        #    others = self.combinations[combElem]
        #    for otherElem in others:
        #        self.addHomologyRelation(combElem, otherElem, {'step': 'finalize', 'relation': (combElem, tuple(others))},)

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

                        #merge props!

                        del self.homologies[xj]
                        changed=True

                        break




    def addCombination(self, id1, listIDs, properties):

        for x in listIDs:
            self.combinations[id1].add(x)

        self.combinationProperties[id1] = mergeDicts(self.combinationProperties.get(id1, None), properties)

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

                if aline[0].startswith('COMBID'):
                    homdb.combinations[aline[0]].add( (aline[1], aline[2]) )
                elif aline[0].startswith('MULCMB'):
                    continue
                elif aline[0].startswith('HOMID'):
                    homdb.homologies[aline[0]].add( (aline[1], aline[2]) )
                else:
                    continue


        return homdb

    def save_to_file(self, saveLocation):

        with open(saveLocation, 'w') as outfile:

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

            outfile.write(outString)

            outStr.close()
            outStr = io.StringIO()

            cnt = 0
            for x in self.combinations:

                elemProps = None
                if x in self.combinationProperties:
                    elemProps = self.combinationProperties[x]

                printPartners = [x] + list(self.combinations[x])

                for partner in printPartners:

                    descriptor = ["COMBID"+str(cnt), str(partner[0]), str(partner[1])]

                    if elemProps != None:
                        descriptor.append( str(elemProps) )

                    outStr.write("\t".join(descriptor) + "\n")

                cnt+= 1

            outfile.write(outStr.getvalue())
            outStr.close()

            outStr = io.StringIO()
            cnt = 0
            for x in self.multiCombinations:
                outStr.write(x.to_db_str("MULCMB" + str(cnt)))
                cnt+=1

            outfile.write(outStr.getvalue())
            outStr.close()

