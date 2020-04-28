#!/usr/bin/env python
import sys
import os
from glob import glob
from cyvcf2 import VCF, Writer
import gzip
import numpy as np

from .sv_phaser import phase_svs
from .snv_phaser import phase_snvs
HOM_ALT= 2
HET = 1
HOM_REF = 0
VCF_TYPES = ["vcf","vcf.gz", "bcf"]
SV_TYPES = ["DEL","DUP","INV", "CNV","DUP:TANDEM","DEL:ME"]
SNV_TYPE = "POINT"

labels = ["chrom","start", "end","kid", "vartype"]

def read_vars_bed(bedname):
    with open(bedname, 'r') as dnms_file:
        for line in dnms_file:
            if line[0] == "#":
                continue
            fields = line.strip().split()
            if not len(fields) == 5:
                sys.exit("dnms bed file must contain the following columns exactly: " + ", ".join(labels))
            vartype = fields[4]
            if vartype not in SV_TYPES:
                vartype = SNV_TYPE

            yield {
                "chrom"     : fields[0],
                "start"     : int(fields[1]),
                "end"       : int(fields[2]),
                "kid"       : fields[3],
                "vartype"   : vartype,
                "bam"       : ""
            }

def read_vars_bedzip(bedzipname):
    with gzip.open(bedzipname, 'rb') as dnms_file:
        for line in dnms_file:
            if line[0] == "#":
                continue
            fields = line.strip().split()
            if not len(fields) == 5:
                sys.exit("dnms bed file must contain the following columns exactly: " + ", ".join(labels))
            
            vartype = fields[4]
            if vartype not in SV_TYPES:
                vartype = SNV_TYPE
            
            yield {
                "chrom"     : fields[0],
                "start"     : int(fields[1]),
                "end"       : int(fields[2]),
                "kid"       : fields[3],
                "vartype"   : vartype,
                "bam"       : ""
            }

def read_vars_vcf(vcfname):
    vcf = VCF(vcfname)
    for variant in vcf:
        chrom = variant.CHROM
        start = variant.start
        end = variant.end
        
        vartype = variant.INFO.get("SVTYPE")
        if vartype is None:
            vartype = SNV_TYPE
        
        for i,gt in enumerate(variant.gt_types):
            if gt in [HET,HOM_ALT]:
                kid = vcf.samples[i]
                yield {
                    "chrom"     : chrom, 
                    "start"     : start,
                    "end"       : end, 
                    "kid"       : kid,
                    "vartype"  : vartype,
                    "bam"       : ""
                }


def get_bam_names(bam_dir, bam_pairs):
    bam_dict = {}
    if bam_dir is not None:
        for bam in glob(os.path.join(bam_dir,"*.bam")):
            sample_id = os.path.splitext(os.path.basename(bam))[0]
            if sample_id not in bam_dict:
                bam_dict[sample_id] = set()
            bam_dict[sample_id].add(bam)
        for cram in glob(os.path.join(bam_dir,"*.cram")):
            sample_id = os.path.splitext(os.path.basename(cram))[0]
            if sample_id not in bam_dict:
                bam_dict[sample_id] = set()
            bam_dict[sample_id].add(cram)
    
    if bam_pairs is not None:
        for bam_pair in bam_pairs:
            sample_id,bam = bam_pair
            if not os.path.exists(bam) or not os.path.isfile(bam):
                sys.exit("invalid filename "+bam)
            #only one match per id using this approach, so overwrite anything entered previously
            bam_dict[sample_id] = set()
            bam_dict[sample_id].add(bam)
    return bam_dict


def parse_ped(ped, kids):
    labels = ["kid","dad", "mom","sex"]
    kid_entries = {}
    with open(ped, 'r') as pedfile:
        for line in pedfile:
            fields = line.strip().split()
            if fields[1] in kids:
                kid_entries[fields[1]] = dict(zip(labels,fields[1:5]))
    return kid_entries

def summarize_record(read_record, include_ambiguous, verbose):
    dad_read_count = len(read_record['dad_readnames'])
    mom_read_count = len(read_record['mom_readnames'])
        
    origin_parent_reads = "NA"
    origin_parent = None
    origin_parent_sites = None
    evidence_count = 0
    other_parent = None
    other_parent_sites = None
    other_parent_reads = "NA"
    evidence_types = "READBACKED"
    
    if (dad_read_count > 0) and (dad_read_count >= 10*mom_read_count ):
        origin_parent = read_record['dad']
        other_parent = read_record['mom']
        evidence_count = dad_read_count
        origin_parent_sites = read_record['dad_sites']
        origin_parent_reads = read_record['dad_reads']
        other_parent_sites = read_record['mom_sites']
        other_parent_reads = read_record['mom_reads']
    elif (mom_read_count > 0) and (mom_read_count >= 10*dad_read_count ):
        origin_parent = read_record['mom']
        other_parent = read_record['dad']
        evidence_count = mom_read_count
        origin_parent_sites = read_record['mom_sites']
        origin_parent_reads = read_record['mom_reads']
        other_parent_sites = read_record['dad_sites']
        other_parent_reads = read_record['dad_reads']
    elif include_ambiguous:
        origin_parent = read_record['dad']+"|"+read_record['mom']
        evidence_count = "{}|{}".format(dad_read_count, mom_read_count)
        origin_parent_sites = read_record['dad_sites']
        origin_parent_reads = read_record['dad_reads']
        other_parent_reads = read_record['mom_reads']
        other_parent_sites = read_record['mom_sites']
        evidence_types = "AMBIGUOUS_READBACKED"
    else:
        return None

    merged_record = {
        'chrom'                 :read_record['region']['chrom'],
        'start'                 :int(read_record['region']['start']),
        'end'                   :int(read_record['region']['end']),
        'vartype'              :read_record['vartype'],
        'kid'                   :read_record['kid'],
        'origin_parent'         :origin_parent,
        'other_parent'          :other_parent,
        'evidence_count'        :evidence_count,
        'evidence_types'        :evidence_types,
    }
    if verbose:
        merged_record['origin_parent_sites'] = origin_parent_sites
        merged_record['origin_parent_reads'] = origin_parent_reads
        merged_record['other_parent_sites'] = other_parent_sites
        merged_record['other_parent_reads'] = other_parent_reads
    return merged_record


def write_vcf_output(in_vcf_name, read_records, include_ambiguous, verbose):
    vcf = VCF(in_vcf_name)
    vcf.add_format_to_header({"ID": "UOP", "Description": "Unfazed-identified origin parent. Paternal:`0`, maternal:`1`, missing:`-1`","Type":'Float', 'Number':'1'})
    vcf.add_format_to_header({"ID": "UOPS", "Description": "Count of pieces of evidence supporing the unfazed-identified origin parent or `-1` if missing","Type":'Float', 'Number':'1'})
    vcf.add_format_to_header({"ID": "UET", "Description": "Unfazed evidence type: `0` (readbacked), `1` (non-readbacked, for CNVs only), 2 (both), or `-1` (missing)","Type":'Float', 'Number':'1'})
    writer = Writer("/dev/stdout", vcf)
    

    sample_dict = dict(zip(vcf.samples, range(len(vcf.samples))))
    for variant in vcf:
        #keys = []
        uop = []
        uops = []
        uet = []
        for i,gt in enumerate(variant.gt_types):
            uop_entry = -1
            uops_entry = -1
            uet_entry = -1

            if gt in [HET, HOM_ALT]:
                vartype = variant.INFO.get("SVTYPE")
                if vartype is None:
                    vartype = SNV_TYPE

                key_fields = {
                    "chrom"     : variant.CHROM,
                    "start"     : variant.start,
                    "end"       : variant.end,
                    "sample"    : vcf.samples[i],
                    "vartype"   : vartype
                }
                key = "{chrom}_{start}_{end}_{sample}_{vartype}".format(**key_fields)
                if key in read_records:
                    record_summary = summarize_record(read_records[key], include_ambiguous, verbose)
                    if record_summary is not None:
                        origin_parent = record_summary['origin_parent']
                        if origin_parent == read_records[key]['dad']:
                            uop_entry = 0
                        elif origin_parent == read_records[key]['mom']:
                            uop_entry = 1

                        uops_entry = record_summary['evidence_count']

                        evidence_types = record_summary['evidence_types'].split(',')
                        uet_entry = -1
                        if "READBACKED" in evidence_types and "NON_READBACKED" in evidence_types:
                            uet_entry = 2
                        elif "READBACKED" in evidence_types:
                            uet_entry = 0
                        elif "NON_READBACKED" in evidence_types:
                            uet_entry = 1

            uop.append(uop_entry)
            uops.append(uops_entry)
            uet.append(uet_entry)

        variant.set_format('UOP',np.array(uop))
        variant.set_format('UOPS',np.array(uops))
        variant.set_format('UET',np.array(uet))

        writer.write_record(variant)


def write_bed_output(read_records, include_ambiguous, verbose):
    header = [
        "#chrom", 
        "start",
        "end",
        "vartype",
        "kid",
        "origin_parent",
        "other_parent",
        "evidence_count",
        "evidence_types"
    ]

    if verbose:
        header += [
            "origin_parent_sites", 
            "origin_parent_reads",
            "other_parent_sites",
            "other_parent_reads"
        ]


    print("\t".join(header))
    template = "\t".join([
        "{chrom}",
        "{start}",
        "{end}",
        "{vartype}",
        "{kid}",
        "{origin_parent}",
        "{other_parent}",
        "{evidence_count}",
        "{evidence_types}"
    ])

    if verbose:
        template += "".join([
        "{origin_parent_sites}",
        "{origin_parent_reads}",
        "{other_parent_sites}",
        "{other_parent_reads}"
    ])
    record_summaries = []

    for key in read_records:
        record_summary = summarize_record(read_records[key], include_ambiguous, verbose)
        if record_summary is not None:
            record_summaries.append(record_summary)
    
    record_summaries = sorted(record_summaries, key = lambda x: (x['chrom'], x['start'], x['end']))
    for mr in record_summaries:
        print(template.format(**mr))


def unfazed(args):
    input_type = ""
    bam_names_dict = get_bam_names(args.bam_dir, args.bam_pairs)
    snvs = []
    svs = []
    reader = ''
    if args.dnms.endswith(".bed"):
        reader = read_vars_bed
        input_type = "bed"
    elif args.dnms.endswith(".bed.gz"):
        reader = read_vars_bedzip
        input_type = "bed"
    elif True in [args.dnms.endswith(vcf_type) for vcf_type in VCF_TYPES]:
        reader = read_vars_vcf
        input_type = "vcf"
    else:
        sys.exit("dnms file type is unrecognized. Must be bed, bed.gz, vcf, vcf.gz, or bcf")


    output_type = args.output_type if args.output_type is not None else input_type
    if output_type == "vcf" and input_type != "vcf":
        print("Invalid option: --output-type is vcf, but input is not a vcf type. "+
                "Rerun with `--output-type bed` or input dnms as one of the following:", 
                ", ".join(VCF_TYPES), 
                file=sys.stderr)
        sys.exit(1)


    kids = []
    for var_fields in reader(args.dnms):
        sample = var_fields["kid"]
        kids.append(sample)
        if not sample in bam_names_dict:
            print("missing alignment file for", sample, file=sys.stderr)
            continue
        elif len(bam_names_dict[sample]) != 1:
            print("multiple alignment files for", sample+".", 
                    "Please specify correct alignment file using --bam-pairs", file=sys.stderr)
        bam = list(bam_names_dict[sample])[0]
        var_fields['bam'] = bam

        if var_fields['vartype'] in SV_TYPES:
            svs.append(var_fields)
        elif var_fields['vartype'] == SNV_TYPE:
            snvs.append(var_fields)
    pedigrees = parse_ped(args.ped, kids)
    phased_svs = {}
    phased_snvs = {}

    if len(svs) > 0:
        phased_svs = phase_svs(svs, kids, pedigrees, args.sites, args.threads)
    if len(snvs) > 0:
        phased_snvs = phase_snvs(snvs, kids, pedigrees,args.sites, args.threads)

    all_phased = {**phased_snvs, **phased_svs}

    if output_type == "vcf":
        write_vcf_output(args.dnms, all_phased, args.include_ambiguous, args.verbose)
    elif output_type == "bed":
        write_bed_output(all_phased, args.include_ambiguous, args.verbose)