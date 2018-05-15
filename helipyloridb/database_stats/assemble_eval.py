from database.genomedb import GenomeDB
from database.homDBAnalyser import HomDBAnalyser
from database.homologydb import HomologyDatabase

if __name__ == '__main__':
    baseDIR = '/mnt/c/dev/data/haas/homdb/'

    genomeDB = GenomeDB(baseDIR + "/genomes", loadAll=False)
    homDB = HomologyDatabase.loadFromFile(baseDIR + "/hpp_comb")

    analyse = HomDBAnalyser(homDB, genomeDB)

    def printHOM(homid):
        print(homid)

        aligned = analyse.cluster_align('HOMID'+str(homid))
        longest = ""
        allseqs = set()

        for rec in sorted(aligned._records, key=lambda x: x.id):

            seq = str(rec.seq).replace('-', '')
            allseqs.add((seq, rec.id))

            if len(seq) > len(longest):
                longest = str(rec.seq).replace('-', '')

            print(rec.seq, rec.id)

        return ('HOMID'+str(homid), longest, set(allseqs))

    allhomids = [1230, 1621, 1418, 1453, 1297, 1841, 1448, 1630, 1649, 1628]
    allhomids = [2027, 1873, 1453, 1418, 1628, 1415, 1439, 1494, 1913, 1423]
    allhomids = [2171, 1694, 1496, 1648, 1902, 1444, 2076, 1297, 1418, 1452]
    allhomids = [1297, 1448, 1230, 1628, 1841, 1630, 1621, 1418, 1649, 1453]
    allhomids = [1448, 1230, 1297, 1649, 1841, 1628, 1418, 1630, 1453, 1621]
    allhomids = [1629, 1848, 1649, 1635, 1448, 1418, 1621, 1459, 1630, 1628]

    allhomids = [1453, 1630, 1621, 1841, 1649, 1448, 1297, 1418, 1628]# w/o 1,15

    allhomids = [1898, 1702, 1649, 1480, 2078, 1484, 1492, 1868, 1452, 1448] + [1907, 1701, 1866, 2026, 1702, 1448, 1492, 1605, 1912, 1902]


    allLongHoms = [1903, 1650, 1692, 1605, 1651, 2027, 2073, 1869, 1863, 1836, 1439, 1910, 1690, 1868, 1866, 1898, 1480, 1452, 2026, 1689, 1621]
    allhomids = allLongHoms

    allhomids = [1621, 2027, 1692]
    allhomids = [1492, 1621, 2027, 1902, 1692]
    allhomids = [1692, 1693, 1338, 1672, 1621, 1492, 1902, 1324, 2027]

    allhomids = [1346, 115, 900, 567, 961, 1534, 670, 851, 608, 96, 402] # allow any
    allhomids = [1345, 1714, 1865, 1377, 1338, 1621, 1935, 1846, 1642, 1612, 1624, 1408, 1361, 1742, 1960, 1414, 1962, 283, 1964, 1684]#tree

    allhomids = [1329, 1668, 1594, 1908, 1673, 1873, 2027]
    allhomids = [1862, 1457, 1515, 1865, 2027, 1876, 1419, 1386, 1455, 1643, 1487, 1673, 2076, 1474, 1773, 1427, 1384, 1434, 1610, 1703, 1668, 1672, 1624, 1742, 1621, 1415, 1649, 1308, 1960, 1213]


    allhomids = [1484, 1773, 1692, 1973, 1324, 1846, 1584, 1875, 2027, 1389, 1858, 1492, 1694, 1303, 1847, 1962, 1643, 1385, 1444, 1960, 1453, 1586, 1564, 1413, 1419, 1479, 1856, 1610, 1338, 1457]
    """
    ['10_N2-085C2', '8_N1-008A2', '4_N1-031C1', '9_N1-045A2', '2_N1-025A2', '7_N1-008A1', '14_1-20A_UB64', '6_Tx30a']
Feature ranking:
1. feature 20 (0.600000) HOMID1484
2. feature 82 (0.400000) HOMID1773
3. feature 107 (0.000000) HOMID1692
4. feature 26 (0.000000) HOMID1973
5. feature 28 (0.000000) HOMID1324
6. feature 29 (0.000000) HOMID1846
7. feature 30 (0.000000) HOMID1584
8. feature 31 (0.000000) HOMID1875
9. feature 32 (0.000000) HOMID2027
10. feature 33 (0.000000) HOMID1389
11. feature 34 (0.000000) HOMID1858
12. feature 35 (0.000000) HOMID1492
13. feature 36 (0.000000) HOMID1694
14. feature 37 (0.000000) HOMID1303
15. feature 38 (0.000000) HOMID1847
16. feature 39 (0.000000) HOMID1962
17. feature 40 (0.000000) HOMID1643
18. feature 41 (0.000000) HOMID1385
19. feature 42 (0.000000) HOMID1444
20. feature 43 (0.000000) HOMID1960
21. feature 44 (0.000000) HOMID1453
22. feature 45 (0.000000) HOMID1586
23. feature 46 (0.000000) HOMID1564
24. feature 47 (0.000000) HOMID1413
25. feature 48 (0.000000) HOMID1419
26. feature 49 (0.000000) HOMID1479
27. feature 50 (0.000000) HOMID1856
28. feature 27 (0.000000) HOMID1610
29. feature 25 (0.000000) HOMID1338
30. feature 52 (0.000000) HOMID1457
HOMID1484 MC 2 NMC 0
HOMID1773 MC 2 NMC 2
HOMID1692 MC 2 NMC 0
HOMID1973 MC 2 NMC 2
HOMID1324 MC 2 NMC 1
HOMID1846 MC 3 NMC 3
HOMID1584 MC 2 NMC 3
HOMID1875 MC 2 NMC 0
HOMID2027 MC 2 NMC 0
HOMID1389 MC 4 NMC 3
HOMID1858 MC 2 NMC 3
HOMID1492 MC 2 NMC 0
HOMID1694 MC 2 NMC 0
HOMID1303 MC 5 NMC 1
HOMID1847 MC 2 NMC 3
HOMID1962 MC 3 NMC 3
HOMID1643 MC 2 NMC 0
HOMID1385 MC 5 NMC 2
HOMID1444 MC 2 NMC 0
HOMID1960 MC 3 NMC 3
HOMID1453 MC 4 NMC 0
HOMID1586 MC 2 NMC 3
HOMID1564 MC 2 NMC 1
HOMID1413 MC 2 NMC 2
HOMID1419 MC 5 NMC 3
HOMID1479 MC 2 NMC 1
HOMID1856 MC 2 NMC 3
HOMID1610 MC 2 NMC 3
HOMID1338 MC 3 NMC 1
HOMID1457 MC 2 NMC 1
[1484, 1773, 1692, 1973, 1324, 1846, 1584, 1875, 2027, 1389, 1858, 1492, 1694, 1303, 1847, 1962, 1643, 1385, 1444, 1960, 1453, 1586, 1564, 1413, 1419, 1479, 1856, 1610, 1338, 1457]
(15, 2)
Training
4_N1-031C1 [1] True True
2_N1-025A2 [1] True True
14_1-20A_UB64 [1] True True
10_N2-085C2 [1] True True
6_Tx30a [0] False True
7_N1-008A1 [0] False True
9_N1-045A2 [0] False True
8_N1-008A2 [0] False True
Prediction
13_N5-004A1 [1] True True
3_N1-029C1 [1] True True
11_N4-029C2 [1] True True
1_N1-024A1 [0] True False
5_N4-016C1 [0] False True
15_N1-024A2 [0] False True
12_N4-050C1 [0] False True
    
    """




    allhomids = [1361, 2027, 1414, 283, 1165]
    #HOMID1361	HOMID1672	HOMID2027	HOMID1414	HOMID283	HOMID1165

    allhomids = [667, 569, 1165, 283, 439, 525]

    allhomids = [1662, 1165, 1241]
    allhomids = [1714, 2027, 1492]

    allhomids = [1692]

    allhomids = [1448, 1453, 1630, 1629, 1621, 1303, 1668, 1515, 1230, 1742, 1635, 1306, 1215] #>=4 and <= 1



    allhomids = [1839, 1795, 1993, 1496, 1429, 1870, 1859, 1603, 1746, 1630, 1324, 1479, 2072, 1794, 2076, 2078, 1668, 1649, 1453, 1692, 1233, 2070, 1459, 1475, 1875, 2027, 1494, 2024, 1696, 2020, 1415, 1694, 1907, 1701, 1491, 2069, 1215, 1297, 1628, 1876, 1700, 1702, 1693, 1484, 1457, 1515, 1188, 1629, 1908, 1874, 2019, 1913, 1688, 1444, 2023, 1267, 1648, 1564, 1766, 1492, 1631, 1676, 1390, 1636, 1621, 1495, 1864, 1672, 1306, 1230, 1848, 1989, 2021, 1478, 1422, 1448, 1984, 1481, 1371, 1902, 1909, 1418, 1742, 1428, 2071, 1748, 1469, 1873, 1635, 1872, 1643, 1841, 1338, 1485, 1705, 1893, 1303, 1487, 1567, 1477]

    #allhomids = [792]

    allhomids = [2027,1621,1692,1693,1672]
    allhomids = [1549, 1753, 1539, 1971, 1585, 1547, 1820, 1545, 1544, 1555, 1787, 1776, 1634]
    allhomids = [1575, 1460, 2031, 1776, 1753, 1750, 1623, 1585, 1556, 933, 1286, 1435, 1558, 1786, 1346, 2081, 1465]

    allhomids = [2105, 1825]
    allhomids = ['933', '1354', '1413', '283', '1792', '2171', '1621']
    allhomids = ['933', '569', '1413', '1621', '1792', '2171', '1354', '1085', '1743']
    allhomids = ['1413', '2090', '2100', '2081', '1743', '1738', '1621', '1792', '1703', '1814', '1338', '1354', '569', '933', '2105', '1085', '1585']

    allhomids = [1515, 1798, 1303, 1306, 1749, 1786, 1230, 1419, 283]


    allhomids = [283, 1792]
    allhomids = [1792, 1823, 1833, 1825, 1821, 1827, 1831, 2172, 1826, 2173, 2176, 2175, 1817, 1820, 1814, 2174, 1835, 2177, 2015, 1818]
    #allhomids = [2017, 1422, 1354, 1749, 2015, 1469, 1675, 1854, 1599, 1786, 1328, 933, 1515, 1747, 1566, 1680, 1677, 1792, 2171, 1165, 1632, 1303, 1772, 1306, 1419, 1798, 1230, 2018, 1794] # 1 == !15

    # new glimmer db
    allhomids = [1915]

    allhomids = [1693, 1192, 1459, 1464, 1454, 1950, 1694, 1683, 1539, 2021, 1714, 1747, 1484, 2030, 1246, 1526, 1968]
    allhomids = [1731, 1694, 1979, 1503, 1526, 1459, 1699, 1968, 1246, 1323, 2030, 1818, 1708, 1760, 1222, 1692, 1235, 1464, 1950, 1683, 1335, 1484, 1539, 1396, 1714, 1693, 1747, 1948, 1535, 1812, 1454, 2021, 1192]


    allhomids = [449]

    allhomids = [int(x) for x in allhomids]



    allseqs = []

    for homid in allhomids:
        allseqs.append(printHOM(homid))

    for elem in allseqs:

        if len(elem[1]) < 50:
            continue


        print(">"+elem[0] + " " + str(len(elem[1])))
        print(elem[1])

        if elem[0] == 'HOMID1286' or elem[0] == 'HOMID2081':
            for idx, seqi in enumerate(elem[2]):
                print(">"+elem[0] + "_" + str(seqi[1]))

                print(seqi[0])
