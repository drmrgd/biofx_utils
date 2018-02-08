#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Protein Domain Retrieval Script
Starting with a correctly formatted HUGO gene ID, retrieve protein domain position
information from EMBL in a JSON format that can be used as a lookup DB in other 
programs. You can either load a comma separated string of IDs, or a batchfile
containing a list of IDs, one per line, to look up.

??Can also output the data as a simple table??
"""
import sys
import os
import argparse
import requests
import json

from pprint import pprint as pp

version = '0.2.121517'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('build_db', 
        help='Build a protein mapping database from EMBL for future annotation.'
            'Will output to a JSON file containing the date.')
    parser.add_argument('-g', '--gene', metavar='<gene1,gene2,gene3...>', 
        help='Comma separated list of gene(s) to look up.')
    parser.add_argument('-f', '--file', metavar='<batchfile>', 
        help='Batchfile of genes, one per line, to look up.')
    parser.add_argument('-o', '--outfile', metavar='<outfile>', 
        help='Write output to file. DEFAULT: STDOUT.')
    parser.add_argument('-v', '--version', action='version', 
        version = '%(prog)s - ' + version )
    args = parser.parse_args()

    gene_list = []
    if args.gene:
        gene_list = args.gene.split(',')
    elif args.file:
        gene_list = proc_batchfile(args.file)
    else:
        sys.stderr.write("ERROR: You must input a list of genes to look up or "
            "a batchfile of genes to look up.\n")
        sys.exit(1)

    if args.outfile:
        outfh = open(outfile, "w")
    else:
        outfh = sys.stdout

    return gene_list, outfh

def proc_batchfile(batchfile):
    with open(batchfile) as fh:
        return [line.rstrip('\n') for line in fh]

def map_uniprot(gene):
    """
    We will get a lot of results if we don't use a specific uniprot accession
    number.  But, juggling that is hard. So, let's have user input normal UCSC
    gene ID, but then map it to uniprot ID before we retrieve results.
    """
    gene_map = {}
    uniprot_db_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'resources/gene_transcript_uniprot_kegg_map.csv'
    )
    try:
        with open(uniprot_db_file) as fh:
            for line in fh:
                elems = line.split(',')
                gene_map[elems[0]] = elems[2]
    except:
        raise
    return gene_map
            
def api_call(url, query):
    session = requests.Session()
    request = session.get(url, params=query)
    try:
        request.raise_for_status()
    except requests.exceptions.HTTPError as error:
        sys.stderr.write(error)
        sys.exit(1)

    json_data = request.json()
    return json_data
    
def main(gene_list, outfh):
    mapped_genes = map_uniprot(gene_list)
    results = {}

    for gene in gene_list:
        try:
            method = mapped_genes[gene]
        except KeyError:
            sys.stderr.write("WARN: Can not find entry for %s in the database "
                "file. Skipping.\n" % gene)
        url = 'https://www.ebi.ac.uk/proteins/api/features/' + method
        query = {
            'categories' : 'DOMAINS_AND_SITES',
        }
        pp(api_call(url, query))
        #results[gene] = api_call(url, query)

    outfile = 'protein_domain_mapping.json'
    sys.stdout.write('Writing results to %s...\n' % outfile)
    sys.stdout.flush()
    with open(outfile, 'w') as outfh:
        json.dump(results, outfh, indent=4)
    sys.stdout.write('Done!\n')


if __name__ == '__main__':
    gene_list, outfh = get_args()
    main(gene_list, outfh)
