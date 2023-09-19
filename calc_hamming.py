#!/usr/bin/env python3
"""
Input a file containing a list of sequences, one per line, and output a CSV file
containing the hamming distance and the levenshtein distance calculation for
each read.
"""
import sys

from pprint import pprint as pp # noqa
from Levenshtein import distance as lev


version = '1.0.20230816'


def calc_hamming(ref_string, next_string):
    if len(ref_string) == len(next_string):
        return sum(el1 != el2 for el1, el2 in zip(ref_string, next_string))


def __splitstr(string):
    return string.split(" ")


def main(counts_file):
    results = {}

    with open(counts_file) as fh:
        # First one is the highest hit rate and will be reference
        ref_string = fh.readline() 
        refseq, refcounts = __splitstr(ref_string.rstrip("\n"))
        results[refseq] = (refcounts, 0, 0)

        for string in fh:
            seq, counts = __splitstr(string.rstrip("\n"))
            results[seq] = (counts, calc_hamming(refseq, seq), lev(refseq, seq))

    with open("edit_distances.csv", "w") as outfh:
        outfh.write(",".join(["sequence", "counts", "hamming_distance", 
                              "levenshtein_distance"]))
        outfh.write("\n")

        for seq, metrics in results.items():
            outmetrics = map(str, metrics)
            outfh.write(",".join([seq, *outmetrics]))
            outfh.write("\n")


if __name__ == "__main__":
    try:
        counts_file = sys.argv[1]
    except IndexError:
        sys.stderr.write("Error: You must input a counts file!")
        sys.exit(1)

    main(counts_file)


