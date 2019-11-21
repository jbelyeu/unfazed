#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
import read_collector
import site_search
MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

def setup_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-b", "--bam", help="BAM/CRAM for kid", default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/crams/1016.cram")
    parser.add_argument("-c", "--chrom", help="regions of DNM SV", default="12", type=str)
    parser.add_argument("-s", "--start", help="regions of DNM SV", default=103653684, type=int)
    parser.add_argument("-e", "--end", help="regions of DNM SV", default=103655220, type=int)

    return parser.parse_args()

def phase(matches):
    parent_of_origin_evidence_counts = {
        'dad': 0,
        'mom': 0
    }

    for ref_alt in matches:
        for match in matches[ref_alt]:
            read = match['read']
            #to avoid issues from indels, use the reference position to index the read
            read_pos = read.get_reference_positions().index(match['pos'])
            if not read_pos:
                continue
            kid_allele = read.query_sequence[read_pos]

            #if the kid base is ref, this read comes from whichever parent has two refs
            #if the kid base is alt, this read comes from whichever parent has two refs
            kid_allele_isref = (kid_allele == match['ref_informative'])
            if (match['dad'].count(kid_allele_isref) == 2):
                parent_of_origin_evidence_counts["dad"] += 1
            else:
                parent_of_origin_evidence_counts["mom"] += 1
    return parent_of_origin_evidence_counts

def main(args):
    region = {
        'chrom' : args.chrom,
        'start' : args.start,
        'end' : args.end,
    }
    reads = read_collector.collect_reads_sv(args.bam, region)

    informative_sites = {
        "1": [],
        "12": [{
            'pos':103653435, 
            'ref_informative': 'A',
            'alt_informative': 'C',
            'kid':[True, False],
            'dad':[True,False], 
            'mom':[True,True]
        }]
    }
    #matches contain all the informative sites that are overlapped by reads, split into ref/alt
    matches = site_search.match_informative_sites(reads, informative_sites, region['chrom'])
    counts = phase(matches)
    print(counts)



if __name__ == "__main__":
    args = setup_args()
    main(args)
