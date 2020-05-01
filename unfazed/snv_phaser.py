#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
from cyvcf2 import VCF
from concurrent.futures import ThreadPoolExecutor,wait


from .read_collector import collect_reads_snv
from .site_searcher import match_informative_sites
from .informative_site_finder import find,get_prefix

MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

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
    region = "{}{}:{}-{}".format(
        get_prefix(vcf_filehandle),
        chrom.strip("chr"),
        pos,
        int(pos)+1)
    for variant in vcf_filehandle(region):
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
    dnm_reads = collect_reads_snv(denovo['bam'], region, denovo['het_sites'], ref,alt)
    matches = match_informative_sites(dnm_reads, informative_sites)

    if len(matches['alt']) <= 0 and len(matches['ref']) <= 0:
        print("No reads overlap informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
        return
    
    counts = phase_by_reads(matches)
    if dad_id in counts:
        dad_informative_sites = [str(c[1]) for c in counts[dad_id]]
        dad_informative_sites = list(set(dad_informative_sites))
        dad_reads = [c[0].query_name for c in counts[dad_id]]
        dad_reads = list(set(dad_reads))
    else:
        dad_informative_sites = []
        dad_reads = []
    
    if mom_id in counts:
        mom_informative_sites = [str(c[1]) for c in counts[mom_id]]
        mom_informative_sites = list(set(mom_informative_sites))
        mom_reads = [c[0].query_name for c in counts[mom_id]]
        mom_reads = list(set(mom_reads))
    else:
        mom_informative_sites = []
        mom_reads = []

    record = {
        'region'            : region,
        'vartype'           : denovo['vartype'],
        'kid'               : denovo['kid'],
        'dad'               : dad_id,
        'mom'               : mom_id,
        'dad_sites'         : dad_informative_sites,
        'mom_sites'         : mom_informative_sites,
        'evidence_type'     : "readbacked",
        'dad_reads'         : dad_reads,
        'mom_reads'         : mom_reads,
        'cnv_dad_sites'     : '',
        'cnv_mom_sites'     : "",
        'cnv_evidence_type' : ""
    }
    key = [str(v) for v in region.values()]+[denovo['kid'], denovo['vartype']]
    records["_".join(key)] = record


def run_read_phasing(dnms, pedigrees, vcf, threads):
    #get informative sites near SNVs for read-backed phasing
    dnms_with_informative_sites = find(dnms, pedigrees, vcf, 5000, threads, whole_region=False)
    records={}
    if threads != 1:
        executor = ThreadPoolExecutor(threads)
        futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        if 'candidate_sites' not in denovo or len(denovo['candidate_sites']) == 0:
            print("No usable informative sites for variant {}:{}-{}".format(
                denovo['chrom'], denovo['start'], denovo['end']), file=sys.stderr)
            continue

        if threads != 1:
            futures.append(executor.submit(
                multithread_read_phasing, 
                denovo, 
                records, 
                vcf, 
                dad_id, 
                mom_id
            ))
        else:
            multithread_read_phasing(denovo, records, vcf, dad_id, mom_id)
    if threads != 1:
        wait(futures)
    return records



#def phase_snvs(args):
def phase_snvs(dnms, kids, pedigrees, sites, threads):
    return run_read_phasing(dnms, pedigrees, sites, threads)

