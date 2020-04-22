#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
import read_collector
import site_searcher
import informative_site_finder
from cyvcf2 import VCF
from concurrent.futures import ThreadPoolExecutor,wait

MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

def setup_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-d", 
        "--dnms", 
        help="bed file of the DNMs of interest,with chrom, start, end, kid_id, bam_location, var_type",
        #default="/Users/jon/Research/scripts/de_novo_sv/unfazed/snvs.bed")
        default="/Users/jon/Research/scripts/de_novo_sv/unfazed/8346.bed")

    parser.add_argument(
        "-v", 
        "--vcf", 
        help="sorted/bgzipped/indexed vcf/bcf file of SNVs to identify informative sites. Must contain each kid and both parents",
        #default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/snv/hg38/18-01-23_WashU-Yandell-CEPH_Sent.Final_1552990128.vcf.gz")
        default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/ceph_37.vcf.gz")

    parser.add_argument(
        "-p", 
        "--ped", 
        help="ped file including the kid and both parent IDs", 
        type=str, 
        default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/16-08-06_WashU-Yandell-CEPH.ped")
    
    parser.add_argument(
        "-t", 
        "--threads", 
        help="number of threads to use", 
        type=int, 
        default=4)


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
    labels = ["chrom","start", "end","kid", "bam","var_type"]
    kids = []
    dnms = []
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            if line[0] == '#': continue
            #TODO add formatting checks to make sure all the necessary fields are present
            dnms.append(dict(zip(labels,line.strip().split()[:6])))
            try:
                dnms[-1][0] = int(dnms[-1][0])
            except KeyError:
                pass
            kids.append(dnms[-1]['kid'])
    return dnms, kids


def phase_by_reads(matches):
    #parent_ids -> list of informative site matches
    origin_parent_data = {}

    for ref_alt in matches:
        for match_info in matches[ref_alt]:
            read = match_info['read']
            for match in match_info['matches']:
                if len(origin_parent_data) == 0:
                    origin_parent_data[match['ref_parent']] = []
                    origin_parent_data[match['alt_parent']] = []
                #to avoid issues from indels, use the reference position to index the read
                try:
                    read_pos = read.get_reference_positions(full_length=True).index(match['pos'])
                except ValueError:
                    continue
                kid_allele = read.query_sequence[read_pos]

                #if the kid base is ref, this read comes from ref_parent
                #if the kid base is alt, this read comes from alt_parent
                #reads are already assigned to ref or alt relative to the de novo
                #kid base is ref or alt relative to the informative site
                read_origin = ""
                if (kid_allele == match['ref_allele']):
                    read_origin = 'ref_parent'
                elif (kid_allele == match['alt_allele']):
                    read_origin = 'alt_parent'
                else:
                    continue

                #ref_parent means the parent that has the reference allele at the informative site
                #ref haplotype means the read comes from the non-denovo haplotye
                #so if the read comes from the ref_parent and the ref haplotype, de novo is on the alt parent
                if (read_origin == "ref_parent"):
                    if (ref_alt == "ref"):
                        origin_parent_data[match['alt_parent']].append([read, match['pos']])
                    else:
                        origin_parent_data[match['ref_parent']].append([read, match['pos']])
                else:
                    if (ref_alt == "ref"):
                        origin_parent_data[match['ref_parent']].append([read, match['pos']])
                    else:
                        origin_parent_data[match['alt_parent']].append([read, match['pos']])
    return origin_parent_data

def get_refalt(chrom, pos, vcf_filehandle, kid_idx):
    alts = []
    ref = None
    for variant in vcf_filehandle("{}{}:{}-{}".format(
            informative_site_finder.get_prefix(vcf_filehandle),
            chrom.strip("chr"),
            pos,
            int(pos)+1)):
        if ref == None:
            ref = variant.REF
        for alt in variant.ALT:
            alts.append(alt)
    return ref,alts

def multithread_read_phasing(denovo, records, vcf, dad_id, mom_id):
    vcf_filehandle = VCF(vcf)
    sample_dict = dict(zip(vcf_filehandle.samples, range(len(vcf_filehandle.samples))))

    region = {
        "chrom" : denovo['chrom'],
        "start" : denovo['start'],
        "end" : denovo['end'],
    }
    if denovo['kid'] not in sample_dict:
        return
    ref,alts = get_refalt(region['chrom'],region['start'],vcf_filehandle, sample_dict[denovo['kid']])
    if len(alts) < 1:
        print("No usable genotype for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
        return
    elif len(alts) > 1:
        print("Too many genotypes for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
        return
    alt = alts[0]
    informative_sites = denovo['candidate_sites']

    #these are reads that support the ref or alt allele of the de novo variant
    dnm_reads = read_collector.collect_reads_snv(denovo['bam'], region, denovo['het_sites'], ref,alt)
    matches = site_searcher.match_informative_sites(dnm_reads, informative_sites)

    if len(matches['alt']) <= 0 and len(matches['ref']) <= 0:
        print("No reads overlap informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
        return

    counts = phase_by_reads(matches)
    
    #did I put ternary operators and list comprehension in the same lines? What's wrong with me?
    dad_informative_sites = list(set([ str(c[1]) for c in counts[dad_id]])) if dad_id in counts else ["NA"]
    dad_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")" for c in counts[dad_id]]))) if dad_id in counts else "NA"
    dad_readnames = list(set([c[0].query_name for c in counts[dad_id]])) if dad_id in counts else False
    mom_informative_sites = list(set([ str(c[1]) for c in counts[mom_id]])) if mom_id in counts else ["NA"]
    mom_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")"for c in counts[mom_id]]))) if mom_id in counts else "NA"
    mom_readnames = list(set([c[0].query_name for c in counts[mom_id]])) if mom_id in counts else False

    record = {
        'region'            : region,
        'var_type'            : denovo['var_type'],
        'kid'               : denovo['kid'],
        'dad'               : dad_id,
        'mom'               : mom_id,
        'dad_site_count'    : len(dad_informative_sites),
        'mom_site_count'    : len(mom_informative_sites),
        'dad_sites'         : ",".join(dad_informative_sites),
        'mom_sites'         : ",".join(mom_informative_sites),
        'evidence_type'     : "readbacked",
        'dad_reads'         : dad_reads,
        'mom_reads'         : mom_reads,
        'dad_readnames'     : dad_readnames,
        'mom_readnames'     : mom_readnames,
    }
    key = list(region.values())+[denovo['kid'], denovo['var_type']]
    records["_".join(key)] = record

def run_read_phasing(dnms, pedigrees, vcf, threads):
    #get informative sites near SNVs for read-backed phasing
    dnms_with_informative_sites = informative_site_finder.find(dnms, pedigrees, vcf, 5000, args.threads, whole_region=False)
    records={}
    executor = ThreadPoolExecutor(threads)
    futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        if 'candidate_sites' not in denovo or len(denovo['candidate_sites']) == 0:
            print("No usable informative sites for variant {}:{}-{}".format(
                denovo['chrom'], denovo['start'], denovo['end']), file=sys.stderr)
            continue
        futures.append(executor.submit(
            multithread_read_phasing, 
            denovo, 
            records, 
            vcf, 
            dad_id, 
            mom_id
        ))
    wait(futures)
    return records



def main(args):
    dnms,kids = parse_bed(args.dnms)
    pedigrees = parse_ped(args.ped, kids)

    header = [
        "#chrom", 
        "start",
        "end",
        "var_type",
        "kid",
        "origin_parent",
        "other_parent",
        "evidence_count",
        "evidence_types",
        "origin_parent_sites", 
        "origin_parent_reads",
        "other_parent_sites",
        "other_parent_reads",
    ]
    print("\t".join(header))
    template = "\t".join([
        "{chrom}",
        "{start}",
        "{end}",
        "{var_type}",
        "{kid}",
        "{origin_parent}",
        "{other_parent}",
        "{evidence_count}",
        "{evidence_types}",
        "{origin_parent_sites}",
        "{origin_parent_reads}",
        "{other_parent_sites}",
        "{other_parent_reads}",

    ])
    read_records = run_read_phasing(dnms, pedigrees, args.vcf, args.threads)
    merged_records = []

    for key in read_records:

        dad_read_count = len(read_records[key]['dad_readnames'])
        mom_read_count = len(read_records[key]['mom_readnames'])
            
        origin_parent_reads = "NA"
        origin_parent = None
        origin_parent_sites = None
        evidence_count = 0
        other_parent = None
        other_parent_sites = None
        other_parent_reads = "NA"
        evidence_types = "READBACKED"
        
        if (dad_read_count > 0) and (dad_read_count >= 10*mom_read_count ):
            origin_parent = read_records[key]['dad']
            other_parent = read_records[key]['mom']
            evidence_count = dad_read_count
            origin_parent_sites = read_records[key]['dad_sites']
            origin_parent_reads = read_records[key]['dad_reads']
            other_parent_sites = read_records[key]['mom_sites']
            other_parent_reads = read_records[key]['mom_reads']
        elif (mom_read_count > 0) and (mom_read_count >= 10*dad_read_count ):
            origin_parent = read_records[key]['mom']
            other_parent = read_records[key]['dad']
            evidence_count = mom_read_count
            origin_parent_sites = read_records[key]['mom_sites']
            origin_parent_reads = read_records[key]['mom_reads']
            other_parent_sites = read_records[key]['dad_sites']
            other_parent_reads = read_records[key]['dad_reads']
        else:
            origin_parent = read_records[key]['dad']+"|"+read_records[key]['mom']
            evidence_count = "{}|{}".format(dad_read_count, mom_read_count)
            origin_parent_sites = read_records[key]['dad_sites']
            origin_parent_reads = read_records[key]['dad_reads']
            other_parent_reads = read_records[key]['mom_reads']
            other_parent_sites = read_records[key]['mom_sites']
            evidence_types = "AMBIGUOUS_READBACKED"

        merged_record = {
            'chrom'                 :read_records[key]['region']['chrom'],
            'start'                 :int(read_records[key]['region']['start']),
            'end'                   :int(read_records[key]['region']['end']),
            'var_type'                :read_records[key]['var_type'],
            'kid'                   :read_records[key]['kid'],
            'origin_parent'         :origin_parent,
            'other_parent'          :other_parent,
            'evidence_count'        :evidence_count,
            'evidence_types'        :evidence_types,
            'origin_parent_sites'   :origin_parent_sites,
            'origin_parent_reads'   :origin_parent_reads,
            'other_parent_sites'   :other_parent_sites,
            'other_parent_reads'   :other_parent_reads,
        }
        merged_records.append(merged_record)
    merged_records = sorted(merged_records, key = lambda x: (x['chrom'], x['start'], x['end']))
    for mr in merged_records:
        print(template.format(**mr))

    

if __name__ == "__main__":
    args = setup_args()
    main(args)
