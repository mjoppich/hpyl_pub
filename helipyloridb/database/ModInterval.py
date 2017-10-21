from intervaltree import Interval


class ModInterval(Interval):

    def __len__(self):
        return self.end-self.begin+1

    def __eq__(self, other):
        if other == None:
            return False

        return super(ModInterval, self).__eq__(other)

    def __hash__(self):
        """
        Depends on begin and end only.
        :return: hash
        :rtype: Number
        """
        return hash((self.begin, self.end))

    def intersection(self, other):

        if not self.overlaps(other):
            return None

        if self < other:
            x1 = self
            x2 = other
        else:
            x1 = other
            x2 = self

        return ModInterval(x2.begin, x1.end)