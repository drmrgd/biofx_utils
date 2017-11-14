#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Input a file, comma separated list, or a single ClinVar ID(s) and get a table of 
variant information derived from ClinVar using the eutils API functionality of
NCBI. No filtering possible for now, but will be added later.
"""
import sys
import os
import json
import requests
import re
import argparse
import csv

from pprint import pprint as pp

version = '2.0_111417'


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position=100, width=150),
        description = __doc__,
    )
    parser.add_argument('clinvar_id', metavar = '<clinvar_id>', nargs = '?', 
        help='ID or Comma separated list of IDs to search')
    parser.add_argument('-d','--delimiter', metavar = '<delimiter>', default = 'tab',
        help='Delimiter to use for output. Can choose from "comma" or "tab" for now') 
    parser.add_argument('-b', '--batch', metavar = '<batchfile>', 
        help='Batchfile of IDs to search')
    parser.add_argument('-o', '--output', metavar = '<output_file>', 
        help = "Output file to write to. DEFAULT: STDOUT")
    parser.add_argument('-v', '--version', action='version', 
        version = '%(prog)s - ' + version)
    args = parser.parse_args()

    if not args.clinvar_id and not args.batch:
        sys.stderr.write("ERROR: You must either enter a single ClinVar ID to search or a batchfile of IDs.\n")
        sys.stderr.write("USAGE: {} -b <{}> | <{}>\n".format(sys.argv[0],'clinvar_id_list.file', 'clinvar_id')) 
        sys.exit(1)

    if args.delimiter not in ('comma','tab'):
        sys.stderr.write('ERROR: You must choose either "tab" or "comma" when using the delmiter option.\n')
        sys.exit(1)

    return args
                
def clinical_significance_filter(data, input_terms):
    """
    Only output results that match the input clincial significance terms as 
    defined in ClinVar.
    """
    # TODO: Work on this filter.
    terms = {
        'Benign' : '',
        'Likely benign' : '',
        'Uncertain significance' : '',
        'Likely pathogenic' : '',
        'Pathogenic' : '',
        'drug response' : '',
        'association' : '',
        'risk factor' : '',
        'protective' : '',
        'Affects' : '',
        'Conflicting interpretations of pathogenicity' : '',
        'clinical_significance' : '',
        'other' : '',
        'not provided' : '',
    }


def read_varlist(varlist):
    with open(varlist) as fh:
        return [x.rstrip('\n') for x in fh if x]

def get_pmid(varid):
    """
    Input the same Clinvar ID and return PMIDs that cite the clinical evidence 
    for that call. This will return the PubMed link under the "Citations" 
    section of the Clinvar report.
    """
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=clinvar&db=pubmed&id=' + varid + '&format=json'
    pmids = []

    r = requests.get(url)
    ret_data = r.json()

    for i in ret_data['linksets']:
        if 'linksetdbs' in i:
            pmids = [str(link) for link in i['linksetdbs'][0]['links']]
    if pmids:
        # Because Excel absolutely sucks(!!!), delimit these with a semi-colon
        return ';'.join(pmids)
    else:
        return '-'

def get_clinvar_data(varid):
    '''sample url: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=65533&retmode=json'''
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=' + varid + '&retmode=json'
    req = requests.get(url)
    json_data = req.json()
    pmid = get_pmid(varid)
    json_data['result'][varid].update({'pmid' : pmid})
    return json_data

def retrieve_elems(data):
    wanted = ('chr', 'start', 'stop', 'ref', 'alt')
    return [data[x] for x in wanted]

def get_function(aa):
    function = {
        '---'      : 'splicesite',
        'Ter'      : 'nonsense',
        'fs'       : 'frameshift',
    }
    #print("testing {}...".format(aa))
    if aa == '---':
        return function[aa]

    annot = re.search('.*(Ter|fs)$',aa)
    #print(annot.groups())
    try:
        #print('match {}!'.format(annot.group(1)))
        return function[annot.group(1)]
    except AttributeError:
        #print('does not match {}'.format(aa))
        return 'missense'

def get_var_ids(refs,wanted_id):
    ref_list = {elem['db_source'] : elem['db_id'] for elem in refs}
    return ref_list[wanted_id]

def gen_output_handle(arg):
    # TODO: Fix this to use the CSV module so that we can better handle commas in fields.
    if arg:
        print("Writing output to {}".format(arg))
        return open(arg,'w')
    else:
        return sys.stdout

def parse_json(json,varid,outfile,delimiter):
    """
    Parse the JSON file and output only the info we need
    """
    try:
        var_data = json['result'][varid]['variation_set'][0]
    except KeyError:
        sys.stderr.write("ERROR: Can not get data for {}. Skipping.\n".format(varid))
        return

    significance = json['result'][varid]['clinical_significance']['description']
    review_status = json['result'][varid]['clinical_significance']['review_status']
    pmid = json['result'][varid]['pmid']
    results = []
    regex = '(NM.*?)\((.*?)\):(c\..*?)(?: \((.*?)\))?$'

    try:
        dbsnp_id = 'rs' + get_var_ids(var_data['variation_xrefs'],'dbSNP')
    except KeyError:
        #sometimes there is no dbSNP entry, which is OK.  Just fill in value with '-'
        dbsnp_id = '-'

    for entry in var_data['variation_loc']:
        if entry['assembly_name'] == 'GRCh37':
            results = retrieve_elems(entry)
            varname = var_data['variation_name'].replace('&gt;', '>')

            matches = re.search(regex,varname)
            (gene,transcript,cds) = matches.group(2,1,3)
            if matches.group(4):
                aa = matches.group(4)
            else:
                aa = '---'
            function = get_function(aa)

            # results = results + [gene,transcript,cds,aa,function,varid,significance]
            results = results + [gene,transcript,cds,aa,function,varid,dbsnp_id,significance,review_status,pmid]
            results = list(map(lambda x: x if x else '---', results))
    print(delimiter.join(results), file=outfile)

if __name__=='__main__':
    args = get_args()
    if args.batch:
        vlist = read_varlist(args.batch)
    else: 
        vlist = args.clinvar_id.split(',')

    delims = { 'tab' : '\t', 'comma' : ',' } 
    out_fh = gen_output_handle(args.output)
    print(delims[args.delimiter].join(
        ['chr', 'start', 'stop', 'ref', 'alt', 'gene', 'transcript', 'cds', 'aa', 
            'functional', 'clinvar_id', 'dbsnp_id', 'clinical_significance', 
            'review_status', 'PMIDs']
        ), file=out_fh)

    for var in vlist:
        clinvar_json = get_clinvar_data(var)
        parse_json(clinvar_json,var,out_fh,delims[args.delimiter])
