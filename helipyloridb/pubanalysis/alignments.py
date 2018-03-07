from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seq1 = 'TTACGTAAGC'
seq2 = 'TATAAT'

align = pairwise2.align.globalmc(seq1, seq2, 3, -2, lambda x,y: -4*y, lambda x,y: -4*y)

for a in align:
    print(format_alignment(*a))

print(align)