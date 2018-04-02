#!/usr/bin/env python
"""
Using the UCSC LiftOver tool and the Python library established to leverage the
UCSC API, input a comma separated list coordinates, or a simple text file of 
coordinates to map.  The query coordinates must be in the form of "chr#:####", 
and if using an input file, it must have only one coordinate per line. Output 
will be returned in the same fashion, with a comma separated list of original 
coord, new coord. At this time, we can only map from hg18 to hg19, and we 
require an internet connection to use the API.  However, later one, we'll 
develop more functionality.
"""
import sys
import os
import argparse
from pyliftover import LiftOver
from pprint import pprint as pp

version = '0.1.020918'

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('coord', nargs='?', metavar="<chr:pos>",
        help='Coordinate to query.')
    parser.add_argument('-f', '--file', metavar='<batchfile>',
        help='Batch file of coords to process.')
    parser.add_argument('-m', '--mapping', metavar='<query:result>',
        default='hg18:hg19', help='Mapping of original coord to desired coord. '
        'Default: %(default)s.')
    parser.add_argument('-c', '--chain', metavar='<chain_file>', 
        help='Chain file from UCSC for local mapping. Will allow for offline '
        'results.  ****  NOT YET IMPLEMENTED  ****.')
    parser.add_argument('-o', '--outfile', metavar='<output_file>',
        help='File to which to write results. Default STDOUT.')
    parser.add_argument('-v', '--version', action='version', 
        version = '%(prog)s - v' + version)
    args = parser.parse_args()

    coord_list = []
    chainfile = None

    if not any((args.coord,args.file)):
        sys.stderr.write("ERROR: You must input either a list of coords to check"
            " or a batch file containing a set of coords.\n")
        sys.exit(1)
    elif args.coord:
        coord_list = args.coord.split(',')
    else:
        coord_list = proc_batchfile(args.file)

    # Set up proper mapping
    orig_assembly, new_assembly = args.mapping.split(':')
    sys.stderr.write("Mapping from {} to {}.\n".format(orig_assembly, 
        new_assembly))

    #Set up and output filehandle
    if args.outfile:
        sys.stderr.write('Writing output to {}.\n'.format(args.outfile))
        outfh = open(args.outfile, 'w')
    else:
        outfh = sys.stdout

    return coord_list, orig_assembly, new_assembly, chainfile, outfh

def proc_batchfile(batchfile):
    with open(batchfile) as fh:
        return [x.rstrip('\n') for x in fh]

def print_results(results, outfh):
    sys.stderr.write('\n')
    outfh.write(','.join(('OrigChr', 'OrigPos', 'NewChr', 'NewPos')))
    outfh.write('\n')

    for r in results:
        outfh.write(','.join(map(str,r[0:4])))
        outfh.write('\n')

def main(coords, orig_assembly, new_assembly, chainfile, outfh):
    # Create a LiftOver object with desired mapping.
    lo = LiftOver(orig_assembly, new_assembly)

    results = []
    for coord in coords:
        try:
            chrom, pos = coord.split(':')
            # No idea why, but pos needs to be an int instead of a str!
            returnval = lo.convert_coordinate(chrom, int(pos))[0]
            results.append((chrom, pos,) + returnval)
        except:
            # Not sure what kinds of errors we can get.  I think if a locus is
            # deleted, we'll get None as a result (which we'll want to handle),
            # but apart from that, not sure what to expect.  
            sys.stderr.write('Offending coord: %s' % coord)
            raise

    print_results(results, outfh)

if __name__=='__main__':
    coord_list, orig_assembly, new_assembly, chainfile, outfh = get_args()
    try:
        main(coord_list, orig_assembly, new_assembly, chainfile, outfh)
    except KeyboardInterrupt:
        sys.exit(9)

