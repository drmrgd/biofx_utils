#!/usr/bin/env python
# 1/29/2021 - D Sims
################################################################################
"""
Input an Ensembl gene ID or a flat file of Ensembl gene IDs, and return the 
HGNC ID and the official gene symbol using the NCI's db2db service.
"""
import sys
import argparse
import requests

from pprint import pprint as pp # noqa

version = '1.0.020221'

def get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'query',
        nargs='?',
        help='Query string to input.'
    )
    parser.add_argument(
        '-f', '--batchfile',
        metavar="<batchfile>",
        help='File of queries, one per line, to input.'
    )
    parser.add_argument(
        '-o', '--outfile',
        metavar='<output_file>',
        help='Write results to an output file rather than stdout.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'{parser.prog} - v{version}'
    )
    return parser.parse_args()


def api_call(url, query):
        s = requests.Session()
        request = s.get(url, params=query)
        try:
            request.raise_for_status()
        except requests.exceptions.HTTPError as error:
            print(f'{error}')
            raise
        return request.json()

def print_results(result, outfile):
    if outfile:
        sys.stderr.write(f"Writing results to '{outfile}'.\n")
        outfh = open(outfile, 'w')
    else:
        outfh = sys.stdout 

    outfh.write('\t'.join(['Ensembl_ID', 'Gene', 'HNGC_ID']))
    outfh.write('\n')

    for r in result:
        outfh.write('\t'.join(
            [r, result[r].get('gene_symbol'), result[r].get('hgnc_id')]
        ))
        outfh.write('\n')

def __read_batch(f):
    with open(f) as fh:
        #  return ','.join([x.rstrip('\n') for x in fh])
        return [x.rstrip('\n') for x in fh]

def __proc_ret_data(d):
    proc = {}
    for ent in d:
        if ent in ('Input', 'TaxonId'):
            continue

        ensid = d[ent].get('InputValue')

        if len(d[ent]['outputs']) == 0:
            sys.stderr.write(f"No data for {ensid}. Skipping...\n")
            gene_symbols = '-'
            hgnc_ids = '-'
        else:
            try:
                gene_symbols = ','.join(d[ent]['outputs'].get('Gene Symbol', ['-']))
                hgnc_ids     = ','.join(d[ent]['outputs'].get('HGNC ID', ['-']))
            except:
                sys.stderr.write(f"offending record: {d[ent]}\n")
                raise

        proc[ensid] = {
            'gene_symbol' : gene_symbols,
            'hgnc_id'     : hgnc_ids
        }
    return proc

def main(query_list, outfile):
    url ='https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json'

    results = {}
    # Need to chop up the queries to be no more than 250, or else we get an 
    # error that the URL is too long.  So, just split these into 250 query
    # chunks and process one at a time.
    for i in range(0, len(query_list), 250):
        chunk = query_list[i:i+250]
        queries = ','.join(chunk)

        params = {
            'method'      : 'db2db',
            'input'       : 'ensemblgeneid',
            'inputValues' : queries,
            'outputs'     : 'genesymbol,hgncid',
            'taxonId'     : '9606'
        }

        ret_data = api_call(url, params)
        results.update(__proc_ret_data(ret_data))

    print_results(results, outfile)

if __name__ == '__main__':
    args = get_args()

    if not any(x for x in (args.query, args.batchfile)):
        sys.stderr.write("ERROR: You must either input a single query or a "
                "batchfile of queries.")
    elif (args.batchfile):
        queries = __read_batch(args.batchfile)
    else:
        queries = args.query.split(',')

    main(queries, args.outfile)
