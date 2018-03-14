from collections import Counter

allCounts = Counter({('L', 'L'): 44016, ('K', 'K'): 33290, ('I', 'I'): 27323, ('A', 'A'): 26531, ('E', 'E'): 26224, ('S', 'S'): 25686, ('G', 'G'): 22994, ('V', 'V'): 21780, ('F', 'F'): 20787, ('N', 'N'): 20584, ('D', 'D'): 17560, ('T', 'T'): 15991, ('Y', 'Y'): 13990, ('Q', 'Q'): 13591, ('R', 'R'): 12961, ('P', 'P'): 12902, ('M', 'M'): 8917, ('H', 'H'): 8047, ('C', 'C'): 4142, ('W', 'W'): 2672, ('I', 'V'): 583, ('V', 'I'): 572, ('T', 'A'): 408, ('A', 'T'): 405, ('A', 'V'): 329, ('S', 'N'): 322, ('V', 'A'): 306, ('N', 'S'): 301, ('N', 'D'): 293, ('D', 'N'): 290, ('R', 'K'): 274, ('K', 'R'): 267, ('K', 'E'): 242, ('E', 'K'): 210, ('K', 'Q'): 207, ('Q', 'K'): 179, ('I', 'T'): 176, ('G', 'S'): 160, ('L', 'F'): 157, ('L', 'I'): 155, ('F', 'L'): 146, ('E', 'D'): 144, ('D', 'E'): 134, ('T', 'I'): 132, ('A', 'S'): 127, ('N', 'K'): 125, ('P', 'S'): 122, ('K', 'N'): 119, ('I', 'M'): 119, ('I', 'L'): 118, ('S', 'G'): 117, ('S', 'A'): 116, ('R', 'H'): 114, ('G', 'E'): 110, ('H', 'R'): 109, ('H', 'N'): 109, ('N', 'H'): 107, ('S', 'P'): 104, ('E', 'G'): 101, ('Y', 'H'): 99, ('H', 'Y'): 94, ('M', 'I'): 94, ('T', 'V'): 81, ('D', 'G'): 79, ('E', 'Q'): 76, ('G', 'D'): 74, ('K', '-'): 71, ('L', 'S'): 70, ('V', 'L'): 69, ('S', 'L'): 67, ('E', '-'): 66, ('R', 'Q'): 65, ('M', 'V'): 65, ('V', 'M'): 64, ('F', 'S'): 64, ('Q', 'E'): 64, ('T', 'S'): 63, ('-', 'K'): 63, ('-', 'M'): 61, ('T', 'K'): 60, ('E', 'A'): 59, ('Q', 'R'): 58, ('N', '-'): 58, ('V', 'T'): 57, ('-', 'N'): 53, ('-', 'L'): 51, ('A', 'E'): 51, ('L', 'V'): 50, ('S', 'T'): 49, ('N', 'T'): 48, ('K', 'T'): 47, ('L', 'P'): 47, ('T', 'N'): 46, ('S', '-'): 43, ('I', 'A'): 42, ('A', '-'): 42, ('T', '-'): 42, ('A', 'G'): 42, ('H', 'Q'): 41, ('Q', '-'): 41, ('M', 'T'): 40, ('G', 'N'): 40, ('A', 'I'): 40, ('-', 'S'): 39, ('T', 'M'): 39, ('L', '-'): 39, ('-', 'E'): 39, ('Q', 'H'): 38, ('C', 'R'): 37, ('-', 'I'): 36, ('-', 'T'): 36, ('P', 'L'): 36, ('I', '-'): 35, ('S', 'F'): 35, ('-', 'G'): 34, ('P', '-'): 33, ('L', 'M'): 32, ('N', 'G'): 32, ('G', 'A'): 30, ('-', 'R'): 30, ('F', '-'): 28, ('R', 'C'): 28, ('Q', 'P'): 28, ('H', 'S'): 27, ('K', 'S'): 27, ('G', '-'): 27, ('M', 'L'): 27, ('E', 'N'): 27, ('R', '-'): 27, ('F', 'Y'): 26, ('-', 'Q'): 26, ('Y', 'C'): 26, ('P', 'Q'): 26, ('D', 'A'): 25, ('D', 'S'): 25, ('-', 'D'): 24, ('V', '-'): 24, ('C', 'Y'): 23, ('S', 'D'): 23, ('K', 'A'): 23, ('S', 'R'): 23, ('-', 'F'): 23, ('W', 'L'): 23, ('A', 'D'): 22, ('Y', 'F'): 22, ('K', 'G'): 22, ('-', 'V'): 21, ('A', 'M'): 21, ('K', 'D'): 21, ('D', '-'): 21, ('A', 'P'): 21, ('H', 'D'): 21, ('N', 'E'): 21, ('S', 'H'): 21, ('G', 'R'): 20, ('F', 'I'): 20, ('S', 'K'): 20, ('M', 'A'): 20, ('L', 'A'): 20, ('I', 'F'): 20, ('-', 'A'): 19, ('S', 'I'): 19, ('P', 'T'): 19, ('N', 'A'): 19, ('N', 'R'): 19, ('L', 'W'): 18, ('R', 'G'): 18, ('P', 'A'): 18, ('A', 'L'): 18, ('M', '-'): 17, ('I', 'N'): 17, ('Q', 'N'): 17, ('T', 'E'): 17, ('N', 'I'): 17, ('A', 'K'): 16, ('R', 'E'): 16, ('T', 'L'): 16, ('-', 'P'): 16, ('-', 'Y'): 15, ('E', 'T'): 15, ('S', 'Q'): 15, ('R', 'S'): 15, ('I', 'S'): 15, ('R', 'N'): 15, ('H', 'K'): 15, ('A', 'Q'): 15, ('Q', 'S'): 15, ('C', 'G'): 14, ('H', '-'): 14, ('H', 'C'): 13, ('K', 'I'): 13, ('G', 'K'): 13, ('D', 'T'): 13, ('V', 'F'): 13, ('A', 'N'): 12, ('D', 'K'): 12, ('T', 'P'): 12, ('Q', 'T'): 12, ('Y', '-'): 12, ('F', 'C'): 12, ('C', 'H'): 12, ('Y', 'D'): 11, ('-', 'H'): 11, ('Q', 'L'): 11, ('N', 'Q'): 11, ('K', 'L'): 11, ('F', 'V'): 11, ('M', 'K'): 11, ('P', 'K'): 11, ('L', 'Q'): 10, ('E', 'R'): 10, ('Y', 'S'): 10, ('H', 'P'): 10, ('T', 'G'): 10, ('S', 'E'): 10, ('L', 'R'): 9, ('L', 'T'): 9, ('I', 'K'): 9, ('P', 'H'): 9, ('S', 'V'): 9, ('C', 'S'): 9, ('D', 'H'): 9, ('K', 'H'): 8, ('E', 'S'): 8, ('T', 'R'): 8, ('Y', 'N'): 8, ('S', 'Y'): 8, ('T', 'D'): 8, ('L', 'Y'): 8, ('K', 'P'): 8, ('R', 'M'): 8, ('G', 'V'): 8, ('L', 'K'): 7, ('C', 'F'): 7, ('C', '-'): 7, ('F', 'A'): 7, ('S', 'C'): 7, ('N', 'Y'): 7, ('V', 'E'): 7, ('D', 'R'): 7, ('H', 'E'): 6, ('R', 'A'): 6, ('I', 'Y'): 6, ('T', 'H'): 6, ('D', 'Y'): 6, ('Q', 'D'): 6, ('G', 'C'): 6, ('G', 'H'): 6, ('H', 'T'): 6, ('D', 'V'): 6, ('T', 'F'): 6, ('D', 'Q'): 6, ('R', 'Y'): 6, ('W', '-'): 6, ('N', 'V'): 6, ('G', 'T'): 6, ('K', 'V'): 6, ('V', 'S'): 6, ('L', 'N'): 5, ('E', 'I'): 5, ('P', 'D'): 5, ('C', 'L'): 5, ('R', 'T'): 5, ('Y', 'R'): 5, ('L', 'C'): 5, ('Q', 'Y'): 5, ('G', 'W'): 5, ('R', 'P'): 5, ('R', 'L'): 5, ('H', 'F'): 5, ('-', 'W'): 5, ('P', 'I'): 5, ('E', 'V'): 5, ('P', 'E'): 4, ('V', 'P'): 4, ('V', 'Y'): 4, ('V', 'K'): 4, ('P', 'F'): 4, ('F', 'P'): 4, ('Q', 'G'): 4, ('F', 'H'): 4, ('T', 'Q'): 4, ('E', 'H'): 4, ('N', 'C'): 4, ('V', 'G'): 4, ('V', 'N'): 4, ('Y', 'L'): 4, ('Y', 'T'): 4, ('R', 'D'): 4, ('I', 'R'): 4, ('A', 'F'): 3, ('P', 'N'): 3, ('A', 'H'): 3, ('W', 'V'): 3, ('E', 'M'): 3, ('R', 'I'): 3, ('Y', 'A'): 3, ('D', 'P'): 3, ('Q', 'C'): 3, ('Y', 'I'): 3, ('Y', 'G'): 3, ('Q', 'F'): 3, ('V', 'D'): 3, ('I', 'E'): 3, ('F', 'K'): 3, ('M', 'R'): 3, ('A', 'R'): 3, ('W', 'G'): 3, ('M', 'S'): 3, ('Y', 'E'): 3, ('C', 'Q'): 3, ('C', 'N'): 3, ('R', 'V'): 3, ('Y', 'P'): 3, ('P', 'V'): 3, ('C', 'V'): 3, ('I', 'Q'): 3, ('I', 'P'): 3, ('W', 'K'): 3, ('W', 'R'): 3, ('Q', 'I'): 3, ('G', 'L'): 3, ('H', 'G'): 3, ('Q', 'A'): 3, ('N', 'M'): 3, ('F', 'N'): 3, ('E', 'P'): 3, ('N', 'L'): 3, ('P', 'R'): 3, ('I', 'C'): 3, ('L', 'E'): 3, ('N', 'P'): 3, ('K', 'M'): 3, ('H', 'L'): 3, ('G', 'F'): 3, ('Y', 'M'): 3, ('T', 'Y'): 2, ('V', 'Q'): 2, ('Y', 'V'): 2, ('Y', 'K'): 2, ('V', 'W'): 2, ('Y', 'Q'): 2, ('M', 'F'): 2, ('K', 'C'): 2, ('F', 'W'): 2, ('-', 'C'): 2, ('C', 'E'): 2, ('N', 'F'): 2, ('I', 'D'): 2, ('G', 'I'): 2, ('C', 'A'): 2, ('F', 'D'): 2, ('G', 'Q'): 2, ('D', 'I'): 2, ('S', 'M'): 2, ('L', 'G'): 2, ('G', 'P'): 2, ('C', 'I'): 2, ('A', 'Y'): 2, ('T', 'C'): 2, ('E', 'F'): 2, ('Y', 'W'): 2, ('F', 'T'): 2, ('F', 'G'): 2, ('H', 'A'): 2, ('D', 'W'): 2, ('K', 'Y'): 2, ('P', 'G'): 2, ('W', 'Q'): 2, ('E', 'L'): 1, ('F', 'R'): 1, ('D', 'F'): 1, ('N', 'W'): 1, ('M', 'C'): 1, ('P', 'C'): 1, ('W', 'Y'): 1, ('V', 'R'): 1, ('C', 'M'): 1, ('F', 'E'): 1, ('W', 'C'): 1, ('K', 'F'): 1, ('M', 'E'): 1, ('M', 'N'): 1, ('M', 'Y'): 1, ('A', 'W'): 1, ('W', 'S'): 1, ('V', 'H'): 1, ('I', 'H'): 1, ('M', 'P'): 1, ('V', 'C'): 1, ('E', 'C'): 1, ('W', 'T'): 1, ('D', 'L'): 1, ('M', 'H'): 1, ('C', 'T'): 1, ('A', 'C'): 1, ('W', 'F'): 1, ('Q', 'W'): 1, ('F', 'M'): 1, ('I', 'W'): 1, ('Q', 'V'): 1, ('C', 'W'): 1})


summedChanges = 0
for key in allCounts:

    if key[0] == 'W':
        print(key, allCounts[key])

        if key[0] != key[1]:
            summedChanges += allCounts[key]

print(summedChanges)


print()
print()
print()

summedChanges = 0
for key in allCounts:

    if key[1] == 'W':
        print(key, allCounts[key])

        if key[0] != key[1]:
            summedChanges += allCounts[key]

print(summedChanges)