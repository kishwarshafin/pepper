START = 0
STOP = 1


class IntervalTree(object):
    def __init__(self, intervals, min_bucket=32, _extent=None):
        if len(intervals) < min_bucket:
            self.intervals = intervals
            self.left = self.right = None
            return

        left, right = _extent or (min(interval[START] for interval in intervals), max(interval[STOP] for interval in intervals))
        center = (left+right)/2.0

        self.intervals = list()
        lefts = list()
        rights = list()

        # build list of left, right, and overlapping intervals
        for interval in intervals:
            if interval[STOP] < center:
                lefts.append(interval)
            elif interval[START] > center:
                rights.append(interval)
            else:
                self.intervals.append(interval)     # overlapping

        # recursively assign tree objects to left and right intervals
        self.left = lefts and IntervalTree(lefts, min_bucket, (left, center)) or None
        self.right = rights and IntervalTree(rights, min_bucket, (center, right)) or None
        self.center = center

    def find(self, start, stop):
        """
        find all elements between (or overlapping) start and stop
        """
        overlapping = [i for i in self.intervals if i[STOP] >= start
                       and i[START] <= stop]

        if self.left and start <= self.center:
            overlapping += self.left.find(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find(start, stop)

        return overlapping

    def contains_interval_subset(self, interval):
        """
        Test if the query interval is fully within an interval in the interval tree
        :param interval: query interval
        :return: is_subset: boolean value corresponding to the test result
        """
        start_query, stop_query = interval
        matches = self.find(start_query, stop_query)

        is_subset = False
        for match in matches:
            start = match[0]
            stop = match[1]

            if (start <= start_query <= stop) and (start <= stop_query <= stop):
                is_subset = True

        return is_subset

    def __contains__(self, interval):
        start, stop = interval
        matches = self.find(start, stop)

        if len(matches) > 0:
            return True
        else:
            return False
