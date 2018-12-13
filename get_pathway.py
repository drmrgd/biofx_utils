#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get pathway(s) for a given gene, or output a set of genes for a given pathway.
"""
import sys
import os
import argparse
import json
import csv

from pprint import pprint as pp
from collections import defaultdict

version = '1.3.121318'
sys_json = os.path.join(os.path.dirname(__file__), 'resources', 'pathways.json')

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('-g', '--gene', metavar='<gene>', help='Gene or '
        'comma separated list of genes to look up.')
    parser.add_argument('-b', '--batchfile', metavar='<batchfile>', 
        help='File containing a list of genes (one per line) to look up.')
    parser.add_argument('-p', '--pathway', metavar='<pathway>', help='Pathway '
            'to search and return all genes. Since pathways have spaces in the '
            'names, need to surround with quotes when inputting.')
    parser.add_argument('-j', '--json', metavar='<JSON>', default=sys_json, 
        help='JSON file containing gene / pathway mapping info. Default: '
        '%(default)s.')
    parser.add_argument('-o', '--output', metavar="<output_file>", help='Output'
        ' file for data. Default: STDOUT')
    parser.add_argument('-v', '--version', action='version', 
            version = '%(prog)s - v' + version)
    args = parser.parse_args()

    if not any(x for x in (args.gene, args.batchfile, args.pathway)):
        sys.stderr.write('ERROR: You must input either a gene or pathway to '
            'query!\n')
        sys.exit(1)
    return args

def parse_json(jfile):
    with open(jfile) as fh:
        return json.load(fh)

def get_pathway_by_gene(pathway_data, gene_list):
    """
    Input a gene list and output a set of pathways that correspond to that 
    mapping.
    """
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
    for k,v in data.items():
        # Just output a long string of junk if no results found so that we can
        # try to add something later. Will remove later.
        outdata = ['?????????????']
        if v:
            outdata = v
        csv_out.writerow([k]+outdata)

def read_batchfile(batchfile):
    with open(batchfile) as fh:
        return [x.rstrip('\n') for x in fh]

def main(genes, pathway, jfile, outfile):
    data = parse_json(jfile)

    if genes:
        results = get_pathway_by_gene(data, genes)
    elif pathway:
        results = get_gene_by_pathway(data, pathway)
    
    if outfile:
        sys.stderr.write('Writing output to {}.\n'.format(outfile))
        outfh = open(outfile, 'w')
        csv_out = csv.writer(outfh, lineterminator=os.linesep)
    else:
        csv_out = csv.writer(sys.stdout, lineterminator=os.linesep)

    print_results(results, csv_out)

if __name__ == '__main__':
    args = get_args()
    genes = []
    if args.batchfile:
        genes = read_batchfile(args.batchfile)
    elif args.gene:
        genes = args.gene.split(',')

    try:
        main(genes, args.pathway, args.json, args.output)
    except KeyboardInterrupt:
        sys.exit(9)
