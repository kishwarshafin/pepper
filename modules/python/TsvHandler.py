import csv
from collections import defaultdict
from modules.python.IntervalTree import IntervalTree
# from sys import path
# import time


class TsvHandler:
    def __init__(self, tsv_file_path):
        self.tsv_file_path = tsv_file_path

    def get_bed_intervals(self):
        tsv_file = open(self.tsv_file_path, 'r')
        reader = csv.reader(tsv_file, delimiter='\t')

        intervals = list()
        for line in reader:
            # print(line)
            chromosome_name, start, stop = line[0:3]

            start = int(start)
            stop = int(stop)

            intervals.append([start, stop])

        tsv_file.close()

        return intervals

    def get_bed_intervals_by_chromosome(self, start_offset=0, stop_offset=0, universal_offset=0):
        tsv_file = open(self.tsv_file_path, 'r')
        reader = csv.reader(tsv_file, delimiter='\t')

        intervals_chromosomal = defaultdict(list)
        for line in reader:
            # print(line)
            chromosome_name, start, stop = line[0:3]

            start = int(start) + start_offset + universal_offset
            stop = int(stop) + stop_offset + universal_offset

            intervals_chromosomal[chromosome_name].append([start, stop])

        tsv_file.close()

        return intervals_chromosomal

    def get_subset_of_bed_intervals(self, start, stop, start_offset=0, stop_offset=0, universal_offset=0):
        tsv_file = open(self.tsv_file_path, 'r')
        reader = csv.reader(tsv_file, delimiter='\t')

        intervals = list()
        for line in reader:
            # print(line)
            chromosome_name, bed_start, bed_stop = line[0:3]

            bed_start = int(bed_start) + start_offset + universal_offset
            bed_stop = int(bed_stop) + stop_offset + universal_offset

            # print(type(start))
            # print(type(bed_start))

            if bed_stop >= start and bed_start <= stop:
                intervals.append([bed_start, bed_stop])

        # for interval in intervals:
        #     print(interval)

        tsv_file.close()

        return intervals


    # def get_int
