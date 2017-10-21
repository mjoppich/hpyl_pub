from intervaltree import Interval

from database.ModInterval import ModInterval


class IDObj(ModInterval):

    def __new__(cls, seqid, begin, end, genome):
        modInt = super(IDObj, cls).__new__(cls, int(begin), int(end), None)
        modInt.seqid = seqid
        modInt.genome = genome

        return modInt

    def idtuple(self):
        return (self.genome, self.seqid)

    def __str__(self):

        return "<IDObj(Interval) {genom} {seqid} ({begin}-{end})/>".format(genom=self.genome, seqid=self.seqid, begin=self.begin, end=self.end)

    def __repr__(self):
        return self.__str__()

class DiamondResult:
    def __init__(self):
        #query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, sequence

        self.query = None
        self.subject = None

        self.identity = None
        self.alignment_length = None
        self.mismatches = 0
        self.gap_opens = 0
        self.evalue = 1.0
        self.bit_score = 0

    def __str__(self):
        allelems = [self.query.seqid, self.query.begin, self.query.end, self.subject.seqid, self.subject.begin, self.subject.end, self.identity, self.alignment_length, self.mismatches, self.gap_opens, self.evalue, self.bit_score]
        return "\t".join( [str(x) for x in allelems] )

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.alignment_length

    @classmethod
    def from_line(cls, line, qgenom, sgenom):

        aline = [x.strip() for x in line.split('\t')]

        if len(aline) < 12:
            print("Invalid line: " + line)
            exit(-1)

        ret = DiamondResult()
        query = IDObj(aline[0], aline[6], aline[7], qgenom)
        subj = IDObj(aline[1], aline[8], aline[9], sgenom)

        ret.identity = float(aline[2]) / 100.0
        ret.alignment_length = int(aline[3])
        ret.mismatches = int(aline[4])
        ret.gap_opens = int(aline[5])
        ret.evalue = float(aline[10])
        ret.bit_score = float(aline[11])

        ret.subject = subj
        ret.query = query

        return ret