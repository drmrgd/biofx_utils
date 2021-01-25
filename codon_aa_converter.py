#!/usr/bin/env python3
# Convert between an amino acid (single or three letter code) and a codon. This
# is a combination of several stand alone scripts I have, plus a bit of an
# overhaul of the way I was doing this before to be a little more robust and
# faster.
#
# 1/25/2021 - D Sims
################################################################################
"""
If input a three base sequence (i.e. codon), return an Amino Acid. If enter
Amino Acid (either three letter or one letter), return a list of possibly codons
that could make that amino acid.
"""
import sys
import argparse

from collections import defaultdict
from pprint import pprint as pp # noqa

version = '1.0.012521'

# Globals
single_letter = list('ACDEFGHIKLMNPQRSTVWY*')
three_letter = ('Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln '
        'Arg Ser Thr Val Trp Tyr Ter').split()

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'query', 
        metavar='<codon | amino_acid>',
        help='String to convert. Either input a 3 base codon or an amino acid '
            'as a single or three letter string. If one would like to convert '
            'more than one string, can input as comma separated list.'
    )
    parser.add_argument(
        '-v', '--version',
        action = 'version',
        version = '%(prog)s - v' + version
    )
    return parser.parse_args()

def convert_aa(query):
    """
    If input single letter amino acid, return three letter. Do the opposite if
    input three letter amino acid code.
    """
    res = ''
    if len(query) == 1:
        conv_dict = dict(zip(single_letter, three_letter))
        res = conv_dict.get(query.title(), None)
    elif len(query) == 3:
        conv_dict = dict(zip(three_letter, single_letter))
        res = conv_dict.get(query.title(), None)
    else:
        sys.stderr.write(f"Error: '{query}' does not seem to be an appropriate "
                "amino acid string.\n")
        sys.exit(1)

    if res is None:
        sys.stderr.write(f"Error: Can't convert {query} to AA.\n")
        sys.exit(1)

    return res

def translate(query, direction):
    """
    Translate the query (either a codon or an amino acid) to the opposite, as
    directed by 'direction'.  Direction can be either 'codon' to translate an
    amino acid to a codon, or 'aa' to translate a codon to an amino acid.
    """
    bases = 'TCAG'
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = ('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVV'
        'AAAADDEEGGGG')
    codon_table = dict(zip(codons, amino_acids))

    # Translate a codon to amino acid
    if direction == 'aa':
        # Allow for U's to be input.  But, have to convert to T first, or else
        # the lookup is a bit big.
        query = query.replace('U', 'T')

        res = codon_table.get(query, None)
        if res is None:
            sys.stderr.write(f"Error: Can not translate codon '{query}' to AA!\n")
            sys.exit(1)

        return (res, convert_aa(res))

    # Translate an amino acid to codon(s).
    if direction == 'codon':
        lookup = defaultdict(list)
        for codon, aa in codon_table.items():
            lookup[aa].append(codon)

        if len(query) == 1:
            res = lookup.get(query, None)
        else:
            res = lookup.get(convert_aa(query), None)

        if res is None:
            sys.stderr.write("Error: Can not translate aa '{query}' to codon!\n")
            sys.exit(1)
        return ','.join(res)

def main(input_string):
    """
    Determine if we have a codon or an amino acid string, and return the
    opposite.
    """
    query_list = input_string.split(',')
    results = []

    for query in query_list:
        if all(x in ('A', 'C', 'T', 'G', 'U') for x in query.upper()):
            # We have a codon and want to convert to an amino acid.
            codon = query.upper()
            #  sys.stderr.write(f'Converting codon {codon} to amino acid.\n')
            aa_single, aa_three = translate(codon, direction='aa')
        else:
            # We have an amino acid and want to convert to a codon.
            if len(query) == 1:
                aa_single = query.upper()
                if aa_single in single_letter:
                    #  sys.stderr.write(f"Converting amino acid '{aa_single}' to "
                        #  "codon.\n")
                    codon = translate(aa_single, direction='codon')
                    aa_three = convert_aa(aa_single)
            elif len(query) == 3:
                aa_three = query.title()
                if aa_three in three_letter:
                    #  sys.stderr.write(f"Converting amino acid '{aa_three}' to "
                        #  "codon.\n")
                    codon = translate(aa_three, direction='codon')
                    aa_single = convert_aa(aa_three)
        results.append((aa_single, aa_three, codon))

    # Print the results.
    print('Single\tThree\tCodon(s)')
    for r in results:
        print('\t'.join(r))

if __name__ == '__main__':
    args = get_args()
    main(args.query)
