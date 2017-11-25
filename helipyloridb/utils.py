from collections import Counter, defaultdict
from itertools import chain

fileLocation = '/home/users/joppich/ownCloud/data/hpyloriDB/'
#fileLocation = '/mnt/c/ownCloud/data/hpyloriDB/'

def mergeCounter( counter1, counter2):

    if counter1 == None:
        return counter2

    if counter2 == None:
        return counter1

    mergedCounter = Counter()

    for x in counter1:
        mergedCounter[x] = counter1[x]

    for x in counter2:
        mergedCounter[x] += counter2[x]

    return mergedCounter

def mergeDefaultDict( dict1, dict2):

    if dict1 == None:
        return dict2

    if dict2 == None:
        return dict1

    mergedDefaultDict = defaultdict(dict1.default_factory)

    for x in dict1:
        mergedDefaultDict[x] = dict1[x]

    for x in dict2:
        if x in mergedDefaultDict:
            mergedDefaultDict[x] = mergedDefaultDict[x] + dict2[x]
        else:
            mergedDefaultDict[x] = dict2[x]

    return mergedDefaultDict

def mergeDicts( dict1, dict2, resultType=dict):
    dict3 = resultType()

    if dict1 == None:
        return dict2

    if dict2 == None:
        return dict1

    for k, v in chain(dict1.items(), dict2.items()):

        if k in dict3:

            if not type(v) == type(dict3[k]):
                raise Exception("You try to merge two different objects!")

            if type(v) == list:

                dict3[k] = dict3[k] + v

            elif type(v) == set:

                dict3[k] = dict3[k].union(v)

            elif type(v) == dict:

                dict3[k] = mergeDicts(dict3[k], v)

            elif type(v) == Counter:

                dict3[k] = mergeCounter(dict3[k], v)

            elif type(v) == defaultdict:
                dict3[k] = mergeDefaultDict(dict3[k], v)

            elif type(v) == int or type(v) == float:

                dict3[k] = dict3[k] + v

            else:

                retSet = set()
                retSet.add(v)
                retSet.add(dict3[k])

                if len(retSet) != 1:
                    dict3[k] = tuple(retSet)
                else:
                    dict3[k] = tuple(retSet)[0]
        else:

            dict3[k] = v

    return dict3