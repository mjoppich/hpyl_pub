from lxml import etree

from utils.utils import fileLocation
from xrefs.TextmineDocument import TextMineDocument
from xrefs.tmSentenicer import Sentenicer


class PfamEntry(TextMineDocument):

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

with open(fileLocation + "/PfamFamily.xml") as fin, open(fileLocation + "/PfamFamily.sent", 'w') as fout:

    parser = etree.XMLParser(remove_blank_text=True)

    tree = etree.parse(fin)

    root = tree.getroot()

    for x in tree.xpath('./entries/*'):

        pfamID = x.attrib.get('acc', None)

        pfamChildren = x.getchildren()

        xComments = x.findall('.//field[@name="comment"]')

        pfamComments = [x.text for x in xComments]

        pfamDesc = [x.text for x in x.getchildren() if x.tag=='description']
        pfamName = [x.text for x in x.getchildren() if x.tag=='name']


        entry = PfamEntry(pfamID, pfamName+pfamDesc, pfamComments)

        ret = makeSents.run(entry)

        for x in ret:
            fout.write(str(x) +"\n")

        #print(pfamID, pfamDesc, pfamComments)

