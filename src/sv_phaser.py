#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
import read_collector
import site_searcher
import informative_site_finder

MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

def setup_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-d", 
        "--dnms", 
        help="bed file of the DNMs of interest,with chrom, start, end, kid_id, bam_location",
        default="dnms.bed")

    parser.add_argument(
        "-v", 
        "--vcf", 
        help="vcf file of SNVs to identify informative sites. Must contain each kid and both parents",
        default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/ceph.bcf")

    parser.add_argument(
        "-p", 
        "--ped", 
        help="ped file including the kid and both parent IDs", 
        type=str, 
        default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/16-08-06_WashU-Yandell-CEPH.ped")

    return parser.parse_args()

def parse_ped(ped, kids):
    labels = ["kid","dad", "mom","sex"]
    kid_entries = {}
    with open(ped, 'r') as pedfile:
        for line in pedfile:
            fields = line.strip().split()
            if fields[1] in kids:
                kid_entries[fields[1]] = dict(zip(labels,fields[1:5]))
    return kid_entries

def parse_bed(bed):
    labels = ["chrom","start", "end","kid", "bam"]
    kids = []
    dnms = []
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            #print(line.strip().split())
            if line[0] == '#': continue
            #TODO add formatting checks to make sure all the necessary fields are present
            dnms.append(dict(zip(labels,line.strip().split()[:5])))
            kids.append(dnms[-1]['kid'])
    return dnms, kids


def phase(matches):
    origin_parent_data = {}

    for ref_alt in matches:
        for match in matches[ref_alt]:
            read = match['read']
            #to avoid issues from indels, use the reference position to index the read
            read_pos = read.get_reference_positions().index(match['pos'])
            if not read_pos:
                continue
            kid_allele = read.query_sequence[read_pos]

            #if the kid base is ref, this read comes from ref_parent
            #if the kid base is alt, this read comes from alt_parent
            if (kid_allele == match['ref_allele']):
                if match['ref_parent'] not in origin_parent_data:
                    origin_parent_data[match['ref_parent']] = []
                origin_parent_data[match['ref_parent']].append([read, match['pos']])

            elif (kid_allele == match['alt_allele']):
                if match['alt_parent'] not in origin_parent_data:
                    origin_parent_data[match['alt_parent']] = []
                origin_parent_data[match['alt_parent']].append([read, match['pos']])
            else:
                print("unknown allele for informative site")
    return origin_parent_data



def main(args):
    dnms,kids = parse_bed(args.dnms)
    pedigrees = parse_ped(args.ped, kids)

    dnms_with_informative_sites = informative_site_finder.find(dnms, pedigrees, args.vcf, 10000)
    
    header = ["chrom", "start","end","kid","dad_id","dad_informative_sites","dad_reads","mom_id","mom_informative_sites","mom_reads"]
    print("#"+"\t".join(header))
    for denovo in dnms_with_informative_sites:
        informative_sites = denovo['candidate_sites']

        region = {
            "chrom" : denovo['chrom'],
            "start" : denovo['start'],
            "end" : denovo['end'],
        }

        if len(informative_sites) <= 0:
            print("No usable informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        #these are reads that support the ref or alt allele of the de novo variant
        dnm_reads = read_collector.collect_reads_sv(denovo['bam'], region)
        matches = site_searcher.match_informative_sites(dnm_reads, informative_sites)

        if len(matches['alt']) <= 0 and len(matches['ref']) <= 0:
            print("No reads overlap informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        counts = phase(matches)
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        dad_informative_sites = ",".join(list(set([ str(c[1]) for c in counts[dad_id]]))) if dad_id in counts else "NA"
        dad_reads = ",".join(list(set([c[0].query_name for c in counts[dad_id]]))) if dad_id in counts else "NA"
        mom_informative_sites = ",".join(list(set([ str(c[1]) for c in counts[mom_id]]))) if mom_id in counts else "NA"
        mom_reads = ",".join(list(set([c[0].query_name for c in counts[mom_id]]))) if mom_id in counts else "NA"

        record = list(region.values()) + [
            denovo['kid'], 
            dad_id,
            dad_informative_sites,
            dad_reads,
            mom_id,
            mom_informative_sites,
            mom_reads
        ]
        print("\t".join(record))


if __name__ == "__main__":
    args = setup_args()
    main(args)
