#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys

from .__init__ import __version__
from .unfazed import unfazed


def pair(arg):
    return [x for x in arg.split(":")]


def setup_args():
    parser = argparse.ArgumentParser(
        prog="unfazed", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed version ({})".format(__version__),
        action="version",
        version="%(prog)s " + str(__version__),
    )
    parser.add_argument(
        "-d",
        "--dnms",
        help="valid VCF OR BED file of the DNMs of interest> If BED,"
        + " must contain chrom, start, end, kid_id, var_type columns",
        required=True,
    )

    parser.add_argument(
        "-s",
        "--sites",
        help="sorted/bgzipped/indexed VCF/BCF file of SNVs to identify "
        + "informative sites. Must contain each kid and both parents",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--ped",
        help="ped file including the kid and both parent IDs",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-b",
        "--bam-dir",
        help="directory where bam/cram files (named {sample_id}.bam or "
        + "{sample_id}.cram) "
        + "are stored for offspring. "
        + "If not included, --bam-pairs must be set",
        type=str,
        required=False,
    )

    parser.add_argument(
        "--bam-pairs",
        help="space-delimited list of pairs in the format {sample_id}:{bam_path} "
        + "where {sample_id} matches an offspring id from the dnm file. "
        + "Can be used with --bam-dir arg, must be used in its absence",
        type=pair,
        nargs="*",
        required=False,
    )

    parser.add_argument(
        "-t", "--threads", help="number of threads to use", type=int, default=2
    )

    parser.add_argument(
        "-o",
        "--output-type",
        help="choose output type. If --dnms is not a VCF/BCF,"
        + " output must be to BED format. Defaults to match --dnms input file",
        type=str,
        choices=["vcf", "bed"],
    )

    parser.add_argument(
        "--include-ambiguous",
        help="include ambiguous phasing results",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--verbose",
        help="print verbose output including sites and reads used for phasing. "
        + "Only applies to BED output",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--outfile",
        help="name for output file. Defaults to stdout",
        default="/dev/stdout",
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="reference fasta file (required for crams)",
        required=False,
    )
    parser.add_argument(
        "--build",
        help="human genome build, used to determine sex "
        + "chromosome pseudoautosomal regions. "
        + "If `na` option is chosen, sex chromosomes will not be auto-phased",
        choices=["37", "38", "na"],
        default="38",
        type=str,
    )
    return parser


def main():
    print("\nUNFAZED v{}".format(__version__), file=sys.stderr)
    parser = setup_args()
    args = parser.parse_args()
    if args.bam_dir is None and args.bam_pairs is None:
        print(
            "\nMissing required argument: --bam-dir or --bam-pairs must be set\n",
            file=sys.stderr,
        )
        sys.exit(parser.print_help())
    unfazed(args)


if __name__ == "__main__":
    sys.exit(main() or 0)
