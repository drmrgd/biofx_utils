#!/usr/bin/env python3
import pysam
import argparse

from pprint import pprint as pp # noqa

version = '1.0.20230815'


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'bam', 
        help="BAM File to process"
    )
    parser.add_argument(
        '-t', '--tag',
        help='Tag to check'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version = f'%(prog)s - v{version}'
    )
    return parser.parse_args()

def main(bam, tag):
    bam_data = pysam.AlignmentFile(bam, 'rb')
    cs_tag_data = [read.get_tag(tag) for read in bam_data.fetch()]

    with open("tags.txt", "w") as outfh:
        for tag in cs_tag_data:
            outfh.write(f"{tag}\n")


if __name__ == "__main__":
    args = get_args()
    main(args.bam, args.tag)

