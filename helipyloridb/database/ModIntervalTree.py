from copy import copy

from intervaltree import IntervalTree, Interval

from database.ModInterval import ModInterval


class ModIntervalTree(IntervalTree):

    def __init__(self, intervals=None):
        super(ModIntervalTree, self).__init__(intervals)


    def sum_intervals(self):

        return sum([len( ModInterval(x.begin, x.end) ) for x in self.all_intervals])

    def addi(self, begin, end, data=None):
        """
        Shortcut for add(Interval(begin, end, data)).

        Completes in O(log n) time.
        """
        return self.add(ModInterval(begin, end, data))

    def chop(self, begin, end, datafunc=None, intervalType=ModInterval):
        """
        Like remove_envelop(), but trims back Intervals hanging into
        the chopped area so that nothing overlaps.
        """
        insertions = set()
        begin_hits = [iv for iv in self[begin] if iv.begin < begin]
        end_hits = [iv for iv in self[end] if iv.end > end]

        if datafunc:
            for iv in begin_hits:
                insertions.add(intervalType(iv.begin, begin, datafunc(iv, True)))
            for iv in end_hits:
                insertions.add(intervalType(end, iv.end, datafunc(iv, False)))
        else:
            for iv in begin_hits:
                insertions.add(intervalType(iv.begin, begin, iv.data))
            for iv in end_hits:
                insertions.add(intervalType(end, iv.end, iv.data))

        self.remove_envelop(begin, end)
        self.difference_update(begin_hits)
        self.difference_update(end_hits)
        self.update(insertions)

    def merge_overlaps(self, data_reducer=None, data_initializer=None, newinttype=ModInterval):
        """
        Finds all intervals with overlapping ranges and merges them
        into a single interval. If provided, uses data_reducer and
        data_initializer with similar semantics to Python's built-in
        reduce(reducer_func[, initializer]), as follows:

        If data_reducer is set to a function, combines the data
        fields of the Intervals with
            current_reduced_data = data_reducer(current_reduced_data, new_data)
        If data_reducer is None, the merged Interval's data
        field will be set to None, ignoring all the data fields
        of the merged Intervals.

        On encountering the first Interval to merge, if
        data_initializer is None (default), uses the first
        Interval's data field as the first value for
        current_reduced_data. If data_initializer is not None,
        current_reduced_data is set to a shallow copy of
        data_initiazer created with
            copy.copy(data_initializer).

        Completes in O(n*logn).
        """
        if not self:
            return

        sorted_intervals = sorted(self.all_intervals)  # get sorted intervals
        merged = []
        # use mutable object to allow new_series() to modify it
        current_reduced = [None]
        higher = None  # iterating variable, which new_series() needs access to

        def new_series():
            if data_initializer is None:
                current_reduced[0] = higher.data
                merged.append(higher)
                return
            else:  # data_initializer is not None
                current_reduced[0] = copy(data_initializer)
                current_reduced[0] = data_reducer(current_reduced[0], higher.data)
                merged.append(newinttype(higher.begin, higher.end, current_reduced[0]))

        for higher in sorted_intervals:
            if merged:  # series already begun
                lower = merged[-1]
                if higher.begin <= lower.end:  # should merge
                    upper_bound = max(lower.end, higher.end)
                    if data_reducer is not None:
                        current_reduced[0] = data_reducer(current_reduced[0], higher.data)
                    else:  # annihilate the data, since we don't know how to merge it
                        current_reduced[0] = None
                    merged[-1] = newinttype(lower.begin, upper_bound, current_reduced[0])
                else:
                    new_series()
            else:  # not merged; is first of Intervals to merge
                new_series()

        self.__init__(merged)