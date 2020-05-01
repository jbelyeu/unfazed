#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor,wait

from .read_collector import collect_reads_sv
from .site_searcher import match_informative_sites
from .informative_site_finder import find

MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

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
    labels = ["chrom","start", "end","kid", "bam","vartype"]
    kids = []
    dnms = []
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            if line[0] == '#': continue
            #TODO add formatting checks to make sure all the necessary fields are present
            dnms.append(dict(zip(labels,line.strip().split()[:6])))
            try:
                dnms[-1][0] = int(dnms[-1][0])
            except:
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
                except:
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


def phase_by_snvs(informative_sites):
    if len(informative_sites) <= 0:
        return None
    #parent_ids -> list of informative site matches
    origin_parent_data = {
        informative_sites[0]['ref_parent'] : [],
        informative_sites[0]['alt_parent'] : []
    }

    #split informative sites up by the parent they say is responsible
    for informative_site in informative_sites:
        origin_parent_data[informative_site[informative_site['kid_allele']]].append(informative_site)
    return origin_parent_data

def multithread_read_phasing(denovo, records, dad_id, mom_id):
    region = {
        "chrom" : denovo['chrom'],
        "start" : denovo['start'],
        "end" : denovo['end'],
    }

    #these are reads that support the ref or alt allele of the de novo variant
    dnm_reads = collect_reads_sv(denovo['bam'], region, denovo['het_sites'])
    matches = match_informative_sites(dnm_reads, denovo['candidate_sites'])

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
    #get informative sites near the breakpoints of SVs for reab-backed phasing
    dnms_with_informative_sites = find(dnms, pedigrees, vcf, 5000, threads, whole_region=False)
    records={}
    if threads != 1:
        executor = ThreadPoolExecutor(threads)
        futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        
        if "candidate_sites" not in denovo:
            continue
        informative_sites = denovo['candidate_sites']

        if 'candidate_sites' not in denovo or len(denovo['candidate_sites']) == 0:
            print("No usable informative sites for variant {}:{}-{}".format(
                denovo['chrom'], denovo['start'], denovo['end']), file=sys.stderr)
            continue
        if threads != 1:
            futures.append(executor.submit(
                multithread_read_phasing, 
                denovo, 
                records, 
                dad_id, 
                mom_id
            ))
        else:
            multithread_read_phasing(denovo, records, dad_id, mom_id)
    if threads != 1:
        wait(futures)
    return records

        
    return records

def multithread_cnv_phasing(denovo, records, dad_id, mom_id):
    region = {
        "chrom" : denovo['chrom'],
        "start" : denovo['start'],
        "end" : denovo['end'],
    }

    origin_data = phase_by_snvs(denovo['candidate_sites'])
    if not origin_data:
        return

    evidence_items = {
        dad_id : [],
        mom_id : []
    }
    for parent in evidence_items:
        if parent in origin_data and len(origin_data[parent]) > 0:
            evidence_items[parent] = [str(o['pos']) for o in origin_data[parent]]
    
    record = {
        'region'            : region,
        'vartype'           : denovo['vartype'],
        'kid'               : denovo['kid'],
        'dad'               : dad_id,
        'mom'               : mom_id,
        'cnv_dad_sites'     : evidence_items[dad_id],
        'cnv_mom_sites'     : evidence_items[mom_id],
        'cnv_evidence_type' : "ALLELE-BALANCE",
        'dad_sites'         : "",
        'mom_sites'         : "",
        'evidence_type'     : "",
        'dad_reads'         : [],
        'mom_reads'         : []
    }
    key = [str(r) for r in region.values()]+[denovo['kid'], denovo['vartype']]
    records["_".join(key)] = record

def run_cnv_phasing(dnms, pedigrees, vcf, threads):
    """
    Specialized phasing for CNVs, using the informative sites from the region with a copy-number change
    """
    #get informative sites inside CNVs for purely SNV-based phasing
    dnms_with_informative_sites = find(dnms, pedigrees, vcf, 0, threads)
    records = {}
    if threads != 1:
        executor = ThreadPoolExecutor(threads)
        futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']

        if denovo['vartype'] not in ['DEL','DUP']:
            continue
    
        if 'candidate_sites' not in denovo or len(denovo['candidate_sites']) == 0:
            print("No usable informative sites for variant {}:{}-{}".format(
                denovo['chrom'], denovo['start'], denovo['end']), file=sys.stderr)
            continue
        if threads != 1:
            futures.append(executor.submit(
                multithread_cnv_phasing, 
                denovo, 
                records, 
                dad_id, 
                mom_id
            ))
        else:
            multithread_cnv_phasing(denovo, records, dad_id, mom_id)
    if threads != 1:
        wait(futures)
    return records


#def phase_svs(args):
def phase_svs(dnms, kids, pedigrees, sites, threads):
    cnv_records = run_cnv_phasing(dnms, pedigrees, sites, threads)
    read_records = run_read_phasing(dnms, pedigrees, sites, threads)
    for key in cnv_records:
        if key not in read_records:
            read_records[key] = cnv_records[key]
        else:
            read_records[key]['cnv_dad_sites'] = cnv_records[key]['cnv_dad_sites']
            read_records[key]['cnv_mom_sites'] = cnv_records[key]['cnv_mom_sites']
            read_records[key]['evidence_type'] += ","+cnv_records[key]['cnv_evidence_type']
    return read_records
