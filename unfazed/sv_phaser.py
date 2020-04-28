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
        
        if "candidate_sites" not in denovo:
            continue
        informative_sites = denovo['candidate_sites']

        if len(informative_sites) <= 0:
            print("No usable informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        #these are reads that support the ref or alt allele of the de novo variant
        dnm_reads = read_collector.collect_reads_sv(denovo['bam'], region, denovo['het_sites'])
        matches = site_searcher.match_informative_sites(dnm_reads, informative_sites)

        if len(matches['alt']) <= 0 and len(matches['ref']) <= 0:
            print("No reads overlap informative sites for variant {chrom}:{start}-{end}".format(**region), file=sys.stderr)
            continue

        counts = phase_by_reads(matches)
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']
        
        #did I put terary operators and list comprehension in the same lines? What's wrong with me?
        dad_informative_sites = list(set([ str(c[1]) for c in counts[dad_id]])) if dad_id in counts else ["NA"]
        dad_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")" for c in counts[dad_id]]))) if dad_id in counts else "NA"
        dad_readnames = list(set([c[0].query_name for c in counts[dad_id]])) if dad_id in counts else False
        mom_informative_sites = list(set([ str(c[1]) for c in counts[mom_id]])) if mom_id in counts else ["NA"]
        mom_reads = ",".join(list(set([c[0].query_name+"("+str(c[0].reference_end)+")"for c in counts[mom_id]]))) if mom_id in counts else "NA"
        mom_readnames = list(set([c[0].query_name for c in counts[mom_id]])) if mom_id in counts else False

        record = {
            'region'            : region,
            'svtype'            : denovo['svtype'],
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

        if 'candidate_sites' not in denovo:
            continue
        origin_data = phase_by_snvs(denovo['candidate_sites'])
        evidence_items = {
            dad_id : ["NA"],
            mom_id : ["NA"]
        }
        if not origin_data:
            continue
        for parent in evidence_items:
            if parent in origin_data and len(origin_data[parent]) > 0:
                evidence_items[parent] = [str(o['pos']) for o in origin_data[parent]]

        record = {
            'region'            : region,
            'svtype'            : denovo['svtype'],
            'kid'               : denovo['kid'],
            'dad'               : dad_id,
            'mom'               : mom_id,
            'dad_site_count'    : len(evidence_items[dad_id]),
            'mom_site_count'    : len(evidence_items[mom_id]),
            'dad_sites'         : ",".join(evidence_items[dad_id]),
            'mom_sites'         : ",".join(evidence_items[mom_id]),
            'evidence_type'     : "NON_READBACKED",
            "origin_parent_reads":  "NA"

        }
        key = [str(r) for r in region.values()]+[denovo['kid'], denovo['svtype']]
        records["_".join(key)] = record
    return records


def main(args):
    dnms,kids = parse_bed(args.dnms)
    pedigrees = parse_ped(args.ped, kids)

    header = [
        "#chrom", 
        "start",
        "end",
        "svtype",
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
        "{svtype}",
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
    cnv_records = run_cnv_phasing(dnms, pedigrees, args.vcf)
    read_records = run_read_phasing(dnms, pedigrees, args.vcf)

    #merge evidences
    merged_records = []
    for key in cnv_records:
        dad_count = cnv_records[key]['dad_site_count']
        mom_count = cnv_records[key]['mom_site_count']
        origin_parent = None
        other_parent = None
        origin_parent_sites = None
        other_parent_sites = None
        evidence_types = "NON_READBACKED"
        if float(dad_count)/mom_count >= .95:
            origin_parent = cnv_records[key]['dad']
            other_parent = cnv_records[key]['mom']
            origin_count = dad_count
            origin_parent_sites = cnv_records[key]['dad_sites']
            other_parent_sites = cnv_records[key]['mom_sites']
        elif float(mom_count)/dad_count >= .95:
            origin_parent = cnv_records[key]['mom']
            other_parent = cnv_records[key]['dad']
            origin_count = mom_count
            origin_parent_sites = cnv_records[key]['mom_sites']
            other_parent_sites = cnv_records[key]['dad_sites']
        else:
            origin_parent = cnv_records[key]['dad']+"|"+cnv_records[key]['mom']
            origin_count = "{}|{}".format(dad_count, mom_count)
            origin_parent_sites = cnv_records[key]['dad_sites']
            other_parent_sites = cnv_records[key]['mom_sites']

        merged_record = {
            'chrom'                 :cnv_records[key]['region']['chrom'],
            'start'                 :int(cnv_records[key]['region']['start']),
            'end'                   :int(cnv_records[key]['region']['end']),
            'svtype'                :cnv_records[key]['svtype'],
            'kid'                   :cnv_records[key]['kid'],
            'origin_parent'         :origin_parent,
            'other_parent'          :other_parent,
            'evidence_count'        :str(origin_count),
            'evidence_types'        :cnv_records[key]['evidence_type'],
            'origin_parent_sites'   :origin_parent_sites,
            'origin_parent_reads'   :"NA",
            'other_parent_sites'    :other_parent_sites,
            'other_parent_reads'    :"NA",
        }

        if key in read_records:
            dad_read_count = len(read_records[key]['dad_readnames'])
            mom_read_count = len(read_records[key]['mom_readnames'])
            
            #dad's the de novo origin
            if (dad_read_count > 0) and (dad_read_count >= 10*mom_read_count ):
                merged_record['origin_parent'] = read_records[key]['dad']
                merged_record['other_parent'] = read_records[key]['mom']
                merged_record['evidence_count'] = str(merged_record['evidence_count'])+","+str(dad_read_count)
                merged_record['origin_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['dad_sites']
                merged_record['origin_parent_reads'] = read_records[key]['dad_reads']
                merged_record['other_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['mom_sites']
                merged_record['other_parent_reads'] = read_records[key]['mom_reads']
                merged_record['evidence_types'] += ",READBACKED"
            #mom's the de novo origin
            elif (mom_read_count > 0) and (mom_read_count >= 10*dad_read_count ):
                merged_record['origin_parent'] = read_records[key]['mom']
                merged_record['other_parent'] = read_records[key]['dad']
                merged_record['evidence_count'] = str(merged_record['evidence_count'])+","+str(mom_read_count)
                merged_record['origin_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['mom_sites']
                merged_record['origin_parent_reads'] = read_records[key]['mom_reads']
                merged_record['other_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['dad_sites']
                merged_record['other_parent_reads'] = read_records[key]['dad_reads']
                merged_record['evidence_types'] += ",READBACKED"
            #phasing failed
            else:
                merged_record['origin_parent'] = read_records[key]['dad']+"|"+read_records[key]['mom']
                merged_record['evidence_count']  += ","+"{}|{}".format(dad_read_count, mom_read_count)
                merged_record['evidence_types'] += ",AMBIGUOUS_READBACKED"
                merged_record['origin_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['dad_sites']
                merged_record['origin_parent_reads'] = read_records[key]['dad_reads']
                merged_record['other_parent_sites'] = merged_record['origin_parent_sites']+","+read_records[key]['mom_sites']
                merged_record['other_parent_reads'] = read_records[key]['mom_reads']
        merged_records.append(merged_record)


    for key in read_records:
        if key in cnv_records:
            continue

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
            'svtype'                :read_records[key]['svtype'],
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
