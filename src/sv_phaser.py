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
    labels = ["chrom","start", "end","kid", "bam","svtype"]
    kids = []
    dnms = []
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            if line[0] == '#': continue
            #TODO add formatting checks to make sure all the necessary fields are present
            dnms.append(dict(zip(labels,line.strip().split()[:6])))
            kids.append(dnms[-1]['kid'])
    return dnms, kids


def phase_by_reads(matches):
    #parent_ids -> list of informative site matches
    origin_parent_data = {}

    for ref_alt in matches:
        for match_info in matches[ref_alt]:
            read = match_info['read']
            for match in match_info['matches']:
                #to avoid issues from indels, use the reference position to index the read
                read_pos = read.get_reference_positions(full_length=True).index(match['pos'])
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


def phase_by_snvs(informative_sites):
    if len(informative_sites) <= 0:
        return None
    #parent_ids -> list of informative site matches
    origin_parent_data = {
        informative_sites[0]['ref_parent'] : [],
        informative_sites[0]['alt_parent'] : []
    }

    #split informative sites up by the parent they say is responsible
    for isite in informative_sites:
        origin_parent_data[isite[isite['kid_allele']]].append(isite)
    return origin_parent_data



def run_read_phasing(dnms, pedigrees, vcf):
    #get informative sites near the breakpoints of SVs for reab-backed phasing
    dnms_with_informative_sites = informative_site_finder.find(dnms, pedigrees, vcf, 5000, whole_region=False)
    records={}
    for denovo in dnms_with_informative_sites:
        region = {
            "chrom" : denovo['chrom'],
            "start" : denovo['start'],
            "end" : denovo['end'],
        }
        
        informative_sites = denovo['candidate_sites']

        if len(informative_sites) <= 0:
            print("No usable informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        #these are reads that support the ref or alt allele of the de novo variant
        dnm_reads = read_collector.collect_reads_sv(denovo['bam'], region)

        matches = site_searcher.match_informative_sites(dnm_reads, informative_sites)

        if len(matches['alt']) <= 0 and len(matches['ref']) <= 0:
            print("No reads overlap informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        counts = phase_by_reads(matches)
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        dad_informative_sites = ",".join(list(set([ str(c[1]) for c in counts[dad_id]]))) if dad_id in counts else "NA"
        dad_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")" for c in counts[dad_id]]))) if dad_id in counts else "NA"
        mom_informative_sites = ",".join(list(set([ str(c[1]) for c in counts[mom_id]]))) if mom_id in counts else "NA"
        mom_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")"for c in counts[mom_id]]))) if mom_id in counts else "NA"
        
        record = {
            'region'            : region,
            'svtype'            : denovo['svtype'],
            'kid'               : denovo['kid'],
            'dad'               : dad_id,
            'mom'               : mom_id,
            'dad_site_count'    : len(dad_informative_sites),
            'mom_site_count'    : len(mom_informative_sites),
            'dad_sites'         : dad_informative_sites,
            'mom_sites'         : mom_informative_sites,
            'evidence_type'     : "readbacked"
        }
        key = list(region.values())+[denovo['kid'], denovo['svtype']]
        records["_".join(key)] = record
    return records


def run_cnv_phasing(dnms, pedigrees, vcf):
    """
    Specialized phasing for CNVs, using the informative sites from the region with a copy-number change
    """
    #get informative sites inside CNVs for purely SNV-based phasing
    dnms_with_informative_sites = informative_site_finder.find(dnms, pedigrees, args.vcf, 0)
    records = {}
    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        if denovo['svtype'] not in ['DEL','DUP']:
            continue

        region = {
            "chrom" : denovo['chrom'],
            "start" : int(denovo['start']),
            "end" : int(denovo['end']),
        }


        origin_data = phase_by_snvs(denovo['candidate_sites'])
        evidence_items = {
            dad_id : "NA",
            mom_id : "NA"
        }
        if not origin_data:
            continue
        for parent in evidence_items:
            if parent in origin_data and len(origin_data[parent]) > 0:
                evidence_items[parent] = ",".join([str(o['pos']) for o in origin_data[parent]])
        
        record = {
            'region'            : region,
            'svtype'            : denovo['svtype'],
            'kid'               : denovo['kid'],
            'dad'               : dad_id,
            'mom'               : mom_id,
            'dad_site_count'    : len(evidence_items[dad_id]),
            'mom_site_count'    : len(evidence_items[mom_id]),
            'dad_sites'         : evidence_items[dad_id],
            'mom_sites'         : evidence_items[mom_id],
            'evidence_type'     : "NON_READBACKED"

        }
        key = [str(r) for r in region.values()]+[denovo['kid'], denovo['svtype']]
        records["_".join(key)] = record
    return records


def main(args):
    dnms,kids = parse_bed(args.dnms)
    pedigrees = parse_ped(args.ped, kids)

    header = ["chrom", "start","end","svtype","kid","origin_parent","evidence_count","evidence_types"]
    print("\t".join(header))
    template = "\t".join(["{chrom}", "{start}","{end}","{svtype}","{kid}","{origin_parent}","{evidence_count}","{evidence_types}"])
    cnv_records = run_cnv_phasing(dnms, pedigrees, args.vcf)
    read_records = run_read_phasing(dnms, pedigrees, args.vcf)

    #merge evidences
    merged_records = []
    for key in cnv_records:
        dad_count = cnv_records[key]['dad_site_count']
        mom_count = cnv_records[key]['mom_site_count']
        origin_parent = None
        if float(dad_count)/mom_count >= .95:
            origin_parent = cnv_records[key]['dad']
            origin_count = dad_count
        elif float(mom_count)/dad_count >= .95:
            origin_parent = cnv_records[key]['mom']
            origin_count = mom_count
        else:
            origin_parent = "AMBIGUOUS"
            origin_count = "{}/{}".format(dad_count, mom_count)

        merged_record = {
            'chrom' :cnv_records[key]['region']['chrom'],
            'start' :int(cnv_records[key]['region']['start']),
            'end'   :int(cnv_records[key]['region']['end']),
            'svtype':cnv_records[key]['svtype'],
            'kid'   :cnv_records[key]['kid'],
            'origin_parent' :origin_parent,
            'evidence_count':origin_count,
            'evidence_types':"NON_READBACKED"
        }

        if key in read_records:
            merged_record['evidence_types'] = merged_record['evidence_types'] + ",READBACKED"
            dad_count = read_records[key]['dad_site_count']
            mom_count = read_records[key]['mom_site_count']
            if 'NA' not in [read_records[key]['dad_sites'],read_records[key]['mom_sites']]:
                merged_record['origin_parent'] = "AMBIGUOUS_readbacked"
            elif dad_count > mom_count:
                merged_record['origin_parent'] = read_records[key]['dad']
                merged_record['evidence_count'] = str(merged_record['evidence_count'])+","+str(dad_count)
            elif mom_count > dad_count:
                merged_record['origin_parent'] = read_records[key]['mom']
                merged_record['evidence_count'] = str(merged_record['evidence_count'])+","+str(mom_count)
        merged_records.append(merged_record)


    for key in read_records:
        if key in cnv_records:
            continue

        dad_count = read_records[key]['dad_site_count']
        mom_count = read_records[key]['mom_site_count']
        origin_parent = None
        evidence_count = 0

        if 'NA' not in [read_records[key]['dad_sites'],read_records[key]['mom_sites']]:
            origin_parent = "AMBIGUOUS_readbacked"
            origin_count = "{}|{}".format(dad_count, mom_count)
        elif dad_count > mom_count:
            origin_parent = read_records[key]['dad']
            evidence_count = dad_count
        elif mom_count > dad_count:
            origin_parent = read_records[key]['mom']
            evidence_count = mom_count

        merged_record = {
            'chrom' :read_records[key]['region']['chrom'],
            'start' :int(read_records[key]['region']['start']),
            'end'   :int(read_records[key]['region']['end']),
            'svtype':read_records[key]['svtype'],
            'kid'   :read_records[key]['kid'],
            'origin_parent' :origin_parent,
            'evidence_count':evidence_count,
            'evidence_types':"READBACKED"
        }
        merged_records.append(merged_record)
    merged_records = sorted(merged_records, key = lambda x: (x['chrom'], x['start'], x['end']))
    for mr in merged_records:
        print(template.format(**mr))

    

if __name__ == "__main__":
    args = setup_args()
    main(args)
