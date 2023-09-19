#!/usr/bin/env python3
"""
Read in a BAM file and get the sequence of soft-clipping in each read.
"""
import re
import sys
import pysam
import argparse

from pprint import pprint as pp # noqa

version = '0.1.080322'


def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'bam', 
        metavar='BAM [FILE]',
        help='BAM File to process'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version = '%(prog)s - v' + version
    )
    return parser.parse_args()

def __get_sc_num(cigar):
    '''
    Read CIGAR and return the number of soft-clipped bases.
    '''
    #  bits = re.findall(r'[0-9]+[DIHNMPS]', cigar)
    #  lens = [x.rstrip("S") for x in bits if "S" in x]
    sc_op = re.search(r'([0-9]+)S$', cigar) # Only last softclipped chunk.
    return sc_op.group(1).rstrip("S") if sc_op is not None else -1

def get_clipped(seq, cigar):
    sc_len = int(__get_sc_num(cigar))
    return seq[len(seq)-sc_len:len(seq)] if sc_len != -1 else ""

def read_bam(bam):
    parsed_bam = pysam.AlignmentFile(bam, "rb")

    reads = []
    for read in parsed_bam:
        collected = {
            'name'   : read.query_name.split('-')[1],
            'mb'     : read.get_tag("MB"),
            'cs'     : read.get_tag("CS"),
            'cigar'  : read.cigarstring,
            'sc_seq' : get_clipped(read.seq, read.cigarstring),
            'seq'    : read.seq
        }
        reads.append(collected)
    return reads

def main(bam):
    reads = read_bam(bam)
    fout = ('name', 'seq', 'sc_seq', 'cigar', 'cs', 'mb')
    print('\t'.join(fout))

    for read in reads:
        out = [read[x] for x in fout]
        print('\t'.join(out))
    
if __name__ == "__main__":
    args = get_args()

    main(args.bam)
