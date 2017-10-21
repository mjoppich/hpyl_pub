from intervaltree import IntervalTree, Interval

tree = IntervalTree()

tree.addi(1, 120)

print(tree)

testInterval = Interval(110, 130)

print(tree.overlaps(testInterval))
print(tree[testInterval])

tree.merge_overlaps()