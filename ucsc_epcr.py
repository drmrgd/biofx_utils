#!/usr/bin/env python
## category Conversion
## desc Converts a list of PCR primer pairs to BED regions
#
# Originally from the bedutils package from the ngsutils package:
#  http://ngsutils.org/modules/bedutils/
# 
# Re-writted ane repurposed for my own ill gains!
#
#  TODO:
#     -  Add an optional name field to the single search query method.
#     -  Fix the TSV input method to accept a primer sequence name
#     -  Try to parallelize in order to speed up batch queries.
#     -  Add in more robust arg parsing?
'''
Converts a list of PCR primer pairs to BED regions

Given a genomic FASTA file and a list of primers,  this will generate a new
FASTA file for the targeted region.  This targeted FASTA file can then be
used for targeted resequencing.

PCR input is expected to be paired FASTA entries with names like:
>primername/1
atatcgtgctacgatc
>primername/2
ttgatcggcatataaa

Tab-delimited input should be:
fwd-primer {tab} rev-primer

Uses http://genome.ucsc.edu/cgi-bin/hgPcr to perform ePCR.

Note: This operates by screen scraping,  so it may break unexpectedly in the
future.

Heavily modified and re-worked by D Sims 082015
'''

import os
import sys
import urllib2
import re
import urllib
import time
from pprint import pprint
from collections import defaultdict


version = '0.8.0_080715'
def usage():
    print __doc__
    print """\
Usage: bedutils fromprimers {opts} -fasta file.fa
       bedutils fromprimers {opts} -tab file.txt
       bedutils fromprimers {opts} fwd_primer rev_primer

Options:
-db         name    DB to use
-perfect    bases   Minimum perfect bases (default: 15)
-good       bases   Minimum good bases (default: 15)
-size       bases   Max product size (default: 4000)
-out        file    File to write the data to (default: stdout)
-flip               Flip the reverse primer (default: False)
-quiet              Repress the detailed output

"""
    sys.exit(1)

def insilico_pcr_tab(fname, db, out, quiet, perfect=15,  good=15,  size=4000,  flip=False ):
    pairs = []
    results = defaultdict(list)

    with open(fname) as fileobj: 
        for line in fileobj:
            fwd, rev = [x.strip().upper() for x in line.strip().split('\t')[0:2]]
            if not (fwd, rev) in pairs:
                pairs.append((fwd, rev))

         # TODO
         # This is a bit broken at the moment.  since we don't require a name as input, there is nothing to pass to the results dict.  can use 'i' for now, but kludgy.
         # I think it would be better to actually add a name to the input like in the FASTA example, built a dict above to pass the data in, and then have a more informative
         # and robust output.
        for i, (fwd, rev) in enumerate(pairs):
            data = insilico_pcr_pairs(fwd, rev, i+1, db, out, quiet, perfect, good, size, flip)
            results[i].append(data)
    return results

def insilico_pcr_fasta(fname, db, out, quiet, perfect=15, good=15, size=4000, flip=False):
    primers = {}
    names = []
    seq = ''
    name = None
    results = defaultdict(list) 

    with open(fname) as fileobj:
        for line in fileobj:
            if line[0] == '>':
                if seq:
                    if not name in primers:
                        primers[name] = [seq, ]
                        names.append(name)
                    else:
                        primers[name].append(seq)

                name = line[1:].strip().split()[0].rsplit('/', 1)[0]
                seq = ''
            else:
                seq += line.strip()

        if not name in primers:
            primers[name] = [seq, ]
        else:
            primers[name].append(seq)

        for name in names:
            data = insilico_pcr_pairs(primers[name][0], primers[name][1], name, db, out, quiet, perfect, good, size, flip )
            results[name].append(data) 
    return results 

def insilico_pcr_pairs(fwd, rev, name, db, out, quiet, perfect=15, good=15, size=4000, flip=False ):
    params = [
        ('db', db),
        ('wp_perfect', perfect),
        ('wp_good', good),
        ('boolshad.wp_flipReverse', '0'),
        ('wp_f', fwd),
        ('wp_r', rev),
        ('Submit', 'submit'),
        ('wp_size', size)
        ]
    if flip:
        params.append(('wp_flipReverse', 'on'))

    attempt = 0
    result = None
    data = []

    while attempt < 5 and not result:
        try:
            #if not quiet:
            sys.stderr.write('Trying %s %s->%s (#%s)...' % (name, fwd, rev, attempt + 1))
            f = urllib2.urlopen('http://genome.ucsc.edu/cgi-bin/hgPcr?%s' % urllib.urlencode(params), timeout=180)
            if f:
                result = f.read()
        except:
            pass

        if not result:
            attempt += 1
            if not quiet:
                sys.stderr.write("Timeout... waiting 30 sec\n")
            time.sleep(30)

    if result:
        count = 0
        if re.search('No matches to ', result): 
            sys.stderr.write('no match for %s ~ %s\n' %(fwd, rev))
            # TODO: need a better way to handle missed results.
            data = ['---','---','---',name, '---', '---', fwd, rev, 'no match!']
        else:
            sys.stderr.write('match!\n')
            for m in re.finditer('><A HREF=".*">(.*):(\d+)([-+])(\d+)</A> (.*?) ([ACTG]+) ([ACTG]+)', result):
                data = [m.group(1), int(m.group(2))-1, m.group(4), name, m.group(3), m.group(5), m.group(6), m.group(7), '']
    return data 

def write_results(data):
    # Write out the data
    if out:
        fh = open(out, 'w')
        sys.stderr.write('Writing to file: %s' % out)
    else:
        fh = sys.stdout

    template = '{0:<8}{4:<18}{1:<9}{2:<13}{3:<13}{5:<9}{6:<9}{7:<35}{8:<35}{9:<}\n'
    fh.write(template.format('line', 'chr', 'start', 'end', 'name', 'strand', 'size', 'forward_primer', 'reverse_primer', ''))
    match_num=0
    for query in data:
        for elem in data[query]:
            fh.write(template.format(match_num+1, *elem)) 
            match_num += 1
    return

if __name__ == '__main__':
    db = None
    perfect = 15
    good = 15
    size = 4000
    flip = False
    fasta = None
    tab = None
    fwd = None
    rev = None
    out = None 
    quiet = False

    last = None
    for arg in sys.argv[1:]:
        if last == '-perfect':
            perfect = int(arg)
            last = None
        elif last == '-good':
            good = int(arg)
            last = None
        elif last == '-size':
            size = int(arg)
            last = None
        elif last == '-db':
            db = arg
            last = None
        elif last == '-fasta' and os.path.exists(arg):
            fasta = arg
            last = None
        elif last == '-tab' and os.path.exists(arg):
            tab = arg
            last = None
        elif last == '-out':
            out = arg
            last = None
        elif arg in ['-out', '-db', '-good', '-perfect', '-size', '-fasta', '-tab']:
            last = arg
        elif arg == '-h':
            usage()
        elif arg == '-flip':
            flip = True
        elif arg == '-quiet':
            quiet = True
        else:
            ok = True
            for base in arg:
                if base not in 'ATGCatgc':
                    ok = False
                    break
            if ok and not fwd:
                fwd = arg
            elif ok:
                rev = arg
            else:
                print "Unknown option: %s" % arg
                usage()
    if not db:
        usage()

    if not fasta and not tab and (not fwd and not rev):
        usage()
    
    if fasta:
        data = insilico_pcr_fasta(fasta, db, out, quiet, perfect, good, size, flip)
    elif tab:
        data = insilico_pcr_tab(tab, db, out, quiet, perfect, good, size, flip)
    else:
        data = defaultdict(list)
        data['stdin_primers'].append(insilico_pcr_pairs(fwd, rev, 'stdin_primers', db, out, quiet, perfect, good, size, flip))

    write_results(data)
