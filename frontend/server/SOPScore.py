from Bio import AlignIO
import Bio.SubsMat.MatrixInfo as submat

class SOPScore:


    def __init__(self, matrix=submat.blosum80):

        self.matrix = matrix

    def get_subst_score(self, aaI, aaJ):

        if (aaI, aaJ) in self.matrix:
            return self.matrix[(aaI, aaJ)]
        elif (aaJ, aaI) in self.matrix:
            return self.matrix[(aaJ, aaI)]

        raise ValueError(str(aaI) +" and " + str(aaJ) + " not in matrix")


    def sop(self, alignment, gap):

        assert( isinstance(alignment, AlignIO.MultipleSeqAlignment) )

        print("Hallo")

        idx = [x for x in range(0, len(alignment))]

        aliLen = alignment.get_alignment_length()
        aliCount = len(alignment)

        allScores = []
        for pos in range(0, aliLen):

            posScore = 0

            for i in range(0, aliCount):
                for j in range(i+1, aliCount):

                    aaI = alignment[i].seq[pos]
                    aaJ = alignment[j].seq[pos]

                    if aaI == '*' and aaJ == '*':
                        continue

                    if aaI == '-' or aaJ == '-':

                        if aaI == aaJ:
                            continue

                        posScore += gap

                    else:

                        try:
                            subScore = self.get_subst_score(aaI, aaJ)

                        except:
                            subScore = -10.0

                        posScore += subScore

            allScores.append(posScore)

        return sum(allScores)






if __name__ == '__main__':

    fin = "./tmps/in.msa"
    alignment = AlignIO.read(fin, "fasta")

    scoring = SOPScore()
    score = scoring.sop(alignment, -5.0)

    print(score)