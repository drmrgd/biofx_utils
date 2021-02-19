#!/usr/bin/env python3
"""
Input a set of BAM files and use samtools to index them all in parallel.
"""
import sys
import os
import shutil
import argparse
import subprocess
import multiprocessing

version = '1.0.021921'
nprocs = 12

def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'bams',
        nargs="+",
        metavar='<BAM files>',
        help="Set of one or more BAM files to be indexed."
    )
    parser.add_argument(
        '-p', '--procs',
        type=int,
        metavar='INT',
        default=nprocs,
        help='Number of processor cores to use. Default: %(default)s.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version = f'%(prog)s - v{version}'
    )
    return parser.parse_args()

def usage():
    return 'USAGE: {} <bam_file(s)>\n'.format(os.path.basename(__file__))
    
def index_bam(bam):
    sys.stdout.write('Indexing {}...\n'.format(bam))
    sys.stdout.flush()
    try:
        subprocess.run(['samtools', 'index', bam], check=True)
    except KeyboardInterrupt:
        sys.exit()
    except:
        sys.stderr.write("  -> ERROR!  Can not index BAM file '{}'. File may "
            "be corrupt\n".format(bam))
        return 

def main(bam_files, procs):
    pool = multiprocessing.Pool(processes=procs)
    try:
        pool.map(index_bam, [bam for bam in bam_files])
    except Exception:
        pool.close()
        pool.join()
        sys.exit(1)
    except KeyboardInterrupt:
        pool.terminate()

if __name__=='__main__':
    args = get_args()

    if shutil.which('samtools') is None:
        sys.stderr.write("ERROR: Can not find the package `samtools` in your "
            "$PATH. Be sure that samtools is installed prior to running.\n")
        sys.exit(1)

    main(args.bams, args.procs)
