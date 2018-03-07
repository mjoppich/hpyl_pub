import os
from shutil import copyfile

import networkx
from jinja2 import Environment, PackageLoader, select_autoescape, Template, FileSystemLoader


class MSAViewer:


    @classmethod
    def jinjaRender(cls,tpl_path, context):
        path, filename = os.path.split(tpl_path)
        return Environment( loader=FileSystemLoader(path or './') ).get_template(filename).render(context)

    @classmethod
    def inputdata2str(cls, data):

        outstr = ""

        for idx, (label, seq) in enumerate(data):

            if idx != 0:
                outstr += "\\n\\\n"

            outstr += ">"+label + "\\n\\\n" + seq

        return outstr

    @classmethod
    def makeMSA(cls, data, location='/tmp/msa.html'):
        """

        :param data: must be a list with tuples (seqname, seq)
        :param location:
        :param name:
        :return:
        """

        this_dir, this_filename = os.path.split(__file__)
        #copyfile(this_dir + '/cytoscape.js', location + "/cytoscape.js")



        mydata = {}
        if isinstance(data, list):
            mydata[0] = cls.inputdata2str(data)
        else:

            for x in data:
                mydata[x] = cls.inputdata2str(data[x])


        content = {
            'allGroups': mydata
        }


        outfile = location
        with open(outfile, 'w') as outfile:

            output = cls.jinjaRender(this_dir + '/msa_viewer_template.html', content)

            outfile.write(output)


        return outfile


if __name__ == '__main__':


    mydata = [("bla1", "ATCG"), ("BLA2", "ATTT")]

    groups = {}
    groups["bla"] = mydata
    groups["blub"] = [("bla3", "ATCG"), ("BLA4", "ATTT")]

    MSAViewer.makeMSA(groups, location="/mnt/c/Users/mjopp/Desktop/msatest/")
