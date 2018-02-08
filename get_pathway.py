#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Get pathway(s) for a given gene, or output a set of genes for a given pathway (
which is a bit tricky).
"""
import sys
import os
import argparse
import json
import csv

from pprint import pprint as pp
from collections import defaultdict

version = '1.0.020818'
sys_json = os.path.join(os.path.dirname(__file__), 'resources', 'pathways.json')

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('-g', '--gene', metavar='<gene>', help='Gene or '
        'comma separated list of genes to look up.')
    parser.add_argument('-p', '--pathway', metavar='<pathway>', help='Pathway '
            'to search and return all genes. Since pathways have spaces in the '
            'names, need to surround with quotes when inputting.')
    parser.add_argument('-j', '--json', metavar='<JSON>', default=sys_json, 
        help='JSON file containing gene / pathway mapping info. Default: '
        '%(default)s.')
    parser.add_argument('-o', '--output', metavar="<output_file>", help='Output'
        ' file for data. Default: STDOUT')
    parser.add_argument('-v', '--version', action='version', version='prog(s) - '
        'v' + version)
    args = parser.parse_args()
    return args

def parse_json(jfile):
    with open(jfile) as fh:
        return json.load(fh)

def get_pathway_by_gene(pathway_data, gene):
    """
    Input a gene list and output a set of pathways that correspond to that mapping.
    """
    gene_list = gene.split(',')
    results = {}
    
    for g in gene_list:
        results[g] = []
        for p in pathway_data.keys():
            if g in pathway_data[p]:
                results[g].append(p)
    return results
        

def get_gene_by_pathway(pathway_data, pathway):
    """
    Input a pathway and return a list of genes indicated in that pathway.
    """
    if pathway == '?':
        pathways = pathway_data.keys()
        sys.stderr.write('Valid pathways are:\n')
        for p in sorted(pathways):
            if p == 'file_info':
                continue
            print('\t%s' % p)
        sys.exit()
    else:
        try:
            return {pathway : pathway_data[pathway]}
        except KeyError:
            sys.stderr.write("ERROR: No such pathway '%s'!\n" % pathway)
            sys.exit(1)

def print_results(data, csv_out):
    for k,v in data.iteritems():
        # Just output a long string of junk if no results found so that we can
        # try to add something later. Will remove later.
        outdata = ['?????????????']
        if v:
            outdata = v
        csv_out.writerow([k]+outdata)

def main(gene, pathway, jfile, outfile):
    data = parse_json(jfile)

    if gene:
        results = get_pathway_by_gene(data, gene)
    elif pathway:
        results = get_gene_by_pathway(data, pathway)
    
    if outfile:
        sys.stderr.write('Writing output to {}.\n'.format(outfile))
        outfh = open(outfile, 'w')
        csv_out = csv.writer(outfh)
    else:
        csv_out = csv.writer(sys.stdout)

    print_results(results, csv_out)

if __name__ == '__main__':
    args = get_args()
    try:
        main(args.gene, args.pathway, args.json, args.output)
    except KeyboardInterrupt:
        sys.exit(9)
