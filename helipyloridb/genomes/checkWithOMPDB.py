import csv

from database.homologydb import HomologyDatabase
from utils import fileLocation

with open(fileLocation + "OMP_Database_v3.csv") as infile:

    ompDB = csv.DictReader(infile)
    homDB = HomologyDatabase.loadFromFile(fileLocation + "hpp12_hp")

    allOrganisms = list(homDB.get_all_organisms())

    testOrgs = set()

    for org in allOrganisms:
        if org in ompDB.fieldnames:
            testOrgs.add(org)

    for line in ompDB:

        ompID = line['OMP_ID']

        if line['OMP_ID'] == 'omp_04':
            line=line

        foundInformation = dict()
        origInformation = dict()

        for org in testOrgs:
            foundName = line[ org ]

            if foundName == None or foundName == '' or foundName == 'NA':
                continue

            origInformation[org] = foundName

        for org in origInformation:

            idtuple = (org, origInformation[org])

            homobj = homDB.findHomologyForID( idtuple )

            if homobj != None:
                foundInformation[org] = homobj

        setFoundHomIDs = set()
        for org in foundInformation:
            setFoundHomIDs.add(foundInformation[org])

        if len(setFoundHomIDs) > 1:
            print("Multiple HOMIDs !", ompID)
            print("ompDB", origInformation)
            print("homID", foundInformation)

        if len(foundInformation) != len(origInformation):
            print("Missing Information in homID!", ompID)
            print("ompDB", origInformation)
            print("homID", foundInformation)

        if len(origInformation) != len(testOrgs):
            """
            this is the case if the omp db does not contain a listing
            """

            newInfo = {}
            for x in origInformation:
                newInfo[x] = origInformation[x]

            for org in testOrgs:
                if org in foundInformation:
                    continue

                setHomIDs = set([ foundInformation[x] for x in foundInformation ])

                newInfoDict = dict()

                orgLocusTags = set()
                for homID in setHomIDs:
                    cluster = homDB.get_homology_cluster( homID )

                    clusterLocTag = cluster.get(org, None)
                    if clusterLocTag != None:
                        orgLocusTags.add( clusterLocTag )
                        newInfoDict[homID] = clusterLocTag

                newInfo[org] = newInfoDict


            #print("More Information Known!", ompID)
            #print("Orig: ", origInformation)
            #print("Found:", newInfo)