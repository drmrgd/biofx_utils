#!/usr/bin/env python
# Quickie to convert HGVSp short to long.  Not fully tested on many edge case
# annotations, but seems to work OK on most simple cases.
#
# 5/12/2022 - D Sims
################################################################################
import sys
import re

from pprint import pprint as pp # noqa

def convert(mut): 
    singles = re.findall(r'([A-Z])', mut)
    triples = [aa_convert(x) for x in singles]
    mapping = dict(zip(singles, triples))

    for aa in mapping.keys():
        pattern = re.compile(r'(?:^p.)?(?!{0}[a-z][a-z]){0}'.format(aa))
        mut = re.sub(pattern, mapping[aa], mut)
    return mut

def aa_convert(query):
    single_letter = list("ACDEFGHIKLMNPQRSTVWY*")
    three_letter = ('Ala Cys Asp Glu Phe Gly His Ile Lys Leu Met Asn Pro Gln '
            'Arg Ser Thr Val Trp Tyr Ter').split()
    mapping = dict(zip(single_letter, three_letter))
    return mapping[query.title()]

def main():
    hgvs_file = sys.argv[1]
    with open(hgvs_file) as fh:
        hgvs_list = [x.rstrip('\n') for x in fh]

    converted = [convert(alt) for alt in hgvs_list]

    for c in converted:
        print(c)

if __name__ == "__main__":
    main()
