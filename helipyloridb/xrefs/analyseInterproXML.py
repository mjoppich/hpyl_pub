from lxml import etree

from utils.utils import fileLocation
from xrefs.TextmineDocument import TextMineDocument
from xrefs.tmSentenicer import Sentenicer


class InterproEntry(TextMineDocument):

    def __init__(self, pfamID, desc, comment):
        super(TextMineDocument, self).__init__()

        self.pfamID = pfamID
        self.desc = desc
        self.comment = comment

        self.allSections = {
            'desc': ". ".join(self.desc),
            'comments': ". ".join(self.comment)
        }

    def abstracts(self):
        return []

    def titles(self):
        return self.desc

    def id(self):
        return self.pfamID

    def name(self):
        return self.pfamID

    def sections(self):
        return self.allSections

    def reflists(self):
        return []

    def texts(self):
        return []

makeSents = Sentenicer()

def getText(elem):
    alltext = [x.strip() for x in elem.itertext() if x.strip() is not ""]

    alltext = [x for x in alltext if len(x) > 1]

    finaltext = " ".join(alltext)

    for x in ['[', ']', '\n']:
        finaltext = finaltext.replace(x, '')

    return finaltext


with open(fileLocation + "/interpro.xml") as fin, open(fileLocation + "/interpro.sent", 'w') as fout:

    parser = etree.XMLParser(remove_blank_text=True)

    tree = etree.parse(fin)

    root = tree.getroot()

    for x in tree.xpath('.//interpro'):

        pfamID = x.attrib.get('id', None)

        pfamChildren = x.getchildren()

        pfamDesc = [x for x in x.getchildren() if x.tag=='abstract']
        pfamName = [x.text for x in x.getchildren() if x.tag=='title']


        pfamDesc = ". ".join([getText(x) for x in pfamDesc] )

        entry = InterproEntry(pfamID, pfamName+[pfamDesc], [])

        ret = makeSents.run(entry)

        for x in ret:
            fout.write(str(x) +"\n")

        #print(pfamID, pfamDesc, pfamComments)

