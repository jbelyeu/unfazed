#!/usr/bin/env python
from __future__ import print_function

import gzip
import os
import sys
from glob import glob

import numpy as np
from cyvcf2 import VCF, Writer

from .snv_phaser import phase_snvs
from .sv_phaser import phase_svs

HOM_ALT = 2
HET = 1
HOM_REF = 0
VCF_TYPES = ["vcf", "vcf.gz", "bcf"]
SV_TYPES = ["DEL", "DUP", "INV", "CNV", "DUP:TANDEM", "DEL:ME", "CPX", "CTX"]
SNV_TYPE = "POINT"

labels = ["chrom", "start", "end", "kid", "vartype"]


def read_vars_bed(bedname):
    with open(bedname, "r") as dnms_file:
        for line in dnms_file:
            if line[0] == "#":
                continue
            fields = line.strip().split()
            if not len(fields) == 5:
                sys.exit(
                    "dnms bed file must contain the following columns exactly: "
                    + ", ".join(labels)
                )
            vartype = fields[4]
            if vartype not in SV_TYPES:
                vartype = SNV_TYPE

            yield {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "kid": fields[3],
                "vartype": vartype,
                "bam": "",
            }


def read_vars_bedzip(bedzipname):
    with gzip.open(bedzipname, "rb") as dnms_file:
        for line in dnms_file:
            if line[0] == "#":
                continue
            fields = line.strip().split()
            if not len(fields) == 5:
                sys.exit(
                    "dnms bed file must contain the following columns exactly: "
                    + ", ".join(labels)
                )

            vartype = fields[4]
            if vartype not in SV_TYPES:
                vartype = SNV_TYPE

            yield {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "kid": fields[3],
                "vartype": vartype,
                "bam": "",
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

        for i, gt in enumerate(variant.gt_types):
            if gt in [HET, HOM_ALT]:
                kid = vcf.samples[i]
                yield {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "kid": kid,
                    "vartype": vartype,
                    "bam": "",
                }


def get_bam_names(bam_dir, bam_pairs, cram_ref):
    bam_dict = {}
    cram_found = False
    if bam_dir is not None:
        for bam in glob(os.path.join(bam_dir, "*.bam")):
            sample_id = os.path.splitext(os.path.basename(bam))[0]
            if sample_id not in bam_dict:
                bam_dict[sample_id] = set()
            bam_dict[sample_id].add(bam)
        for cram in glob(os.path.join(bam_dir, "*.cram")):
            cram_found = True
            sample_id = os.path.splitext(os.path.basename(cram))[0]
            if sample_id not in bam_dict:
                bam_dict[sample_id] = set()
            bam_dict[sample_id].add(cram)

    if bam_pairs is not None:
        for bam_pair in bam_pairs:
            sample_id, bam = bam_pair
            if not os.path.exists(bam) or not os.path.isfile(bam):
                sys.exit("invalid filename " + bam)
            # only one match per id using this approach,#
            # so overwrite anything entered previously
            bam_dict[sample_id] = set()
            bam_dict[sample_id].add(bam)
            if bam[-4:] == "cram":
                cram_found = True

    if cram_found:
        if cram_ref is None:
            sys.exit("Missing reference file for CRAM")
        elif not os.path.isfile(cram_ref):
            sys.exit("Reference file is not valid")
    return bam_dict


def parse_ped(ped, kids):
    labels = ["kid", "dad", "mom", "sex"]
    kid_entries = {}
    missing_parents = []
    with open(ped, "r") as pedfile:
        for line in pedfile:
            fields = line.strip().split()
            if fields[1] in kids:
                if fields[2] == "0" or fields[3] == "0":
                    print(
                        "Parent of sample {} missing from pedigree file, will be skipped".format(
                            fields[1]
                        ),
                        file=sys.stderr,
                    )
                    missing_parents.append(fields[1])
                    continue
                kid_entries[fields[1]] = dict(zip(labels, fields[1:5]))

    for sample in kids:
        if (sample not in kid_entries) and (sample not in missing_parents):
            print(
                "{} missing from pedigree file, will be skipped".format(sample),
                file=sys.stderr,
            )
    return kid_entries


def summarize_autophased(read_record):
    chrom = read_record["region"]["chrom"]
    if chrom.lower().strip("chr") == "y":
        origin_parent = read_record["dad"]
        other_parent = read_record["mom"]
    else:
        origin_parent = read_record["mom"]
        other_parent = read_record["dad"]

    record = {
        "chrom": chrom,
        "start": int(read_record["region"]["start"]),
        "end": int(read_record["region"]["end"]),
        "vartype": read_record["vartype"],
        "kid": read_record["kid"],
        "origin_parent": origin_parent,
        "other_parent": other_parent,
        "evidence_count": 1,
        "evidence_types": ["SEX-CHROM"],
    }
    return record


def summarize_record(read_record, include_ambiguous, verbose):
    if read_record["evidence_type"] == "SEX-CHROM":
        return summarize_autophased(read_record)
    dad_read_count = len(read_record["dad_reads"])
    mom_read_count = len(read_record["mom_reads"])
    origin_parent = None
    origin_parent_sites = []
    origin_parent_reads = []
    evidence_count = 0
    other_parent = None
    other_parent_sites = []
    other_parent_reads = []
    evidence_types = []
    ambig = False

    # logic for readbacked phasing
    if (dad_read_count > 0) and (dad_read_count >= 10 * mom_read_count):
        origin_parent = read_record["dad"]
        other_parent = read_record["mom"]
        evidence_count = len(read_record["dad_sites"])
        origin_parent_sites += read_record["dad_sites"]
        origin_parent_reads += read_record["dad_reads"]
        other_parent_sites += read_record["mom_sites"]
        other_parent_reads += read_record["mom_reads"]
        evidence_types.append("READBACKED")
    elif (mom_read_count > 0) and (mom_read_count >= 10 * dad_read_count):
        origin_parent = read_record["mom"]
        other_parent = read_record["dad"]
        evidence_count = len(read_record["mom_sites"])
        origin_parent_sites += read_record["mom_sites"]
        origin_parent_reads += read_record["mom_reads"]
        other_parent_sites += read_record["dad_sites"]
        other_parent_reads += read_record["dad_reads"]
        evidence_types.append("READBACKED")
    elif dad_read_count > 0 and mom_read_count > 0:
        origin_parent = read_record["dad"] + "|" + read_record["mom"]
        evidence_count = dad_read_count + mom_read_count
        origin_parent_sites += read_record["dad_sites"]
        origin_parent_reads += read_record["dad_reads"]
        other_parent_sites += read_record["mom_sites"]
        other_parent_reads += read_record["mom_reads"]
        evidence_types.append("AMBIGUOUS_READBACKED")
        ambig = True

    # logic for cnv phasing
    dad_cnv_site_count = len(read_record["cnv_dad_sites"])
    mom_cnv_site_count = len(read_record["cnv_mom_sites"])
    if (dad_cnv_site_count > 0) and (dad_cnv_site_count >= 10 * mom_cnv_site_count):
        if origin_parent == read_record["mom"] and not ("READBACKED" in evidence_types):
            # this just became ambiguous because of contradictory results
            origin_parent = None
            evidence_count += dad_cnv_site_count + mom_cnv_site_count
            origin_parent_sites += read_record["cnv_dad_sites"]
            other_parent_sites = read_record["cnv_mom_sites"]
            evidence_types = ["AMBIGUOUS_BOTH"]
            ambig = True
        else:
            # dad is origin
            origin_parent = read_record["dad"]
            other_parent = read_record["mom"]
            evidence_count = dad_cnv_site_count
            origin_parent_sites += read_record["cnv_dad_sites"]
            origin_parent_reads += read_record["dad_reads"]
            other_parent_sites += read_record["mom_sites"]
            other_parent_reads += read_record["mom_reads"]
            if ("AMBIGUOUS_READBACKED" in evidence_types):
                evidence_types.remove("AMBIGUOUS_READBACKED")
            evidence_types.append("ALLELE-BALANCE")

    elif (mom_cnv_site_count > 0) and (mom_cnv_site_count >= 10 * dad_cnv_site_count):
        if (origin_parent == read_record["dad"]) and not ("READBACKED" in evidence_types):
            # this just became ambiguous because of contradictory results
            origin_parent = None
            evidence_count += dad_cnv_site_count + mom_cnv_site_count
            origin_parent_sites += read_record["cnv_dad_sites"]
            other_parent_sites += read_record["cnv_mom_sites"]
            evidence_types = ["AMBIGUOUS_BOTH"]
            ambig = True
        else:
            # mom is origin
            origin_parent = read_record["mom"]
            other_parent = read_record["dad"]
            evidence_count = mom_cnv_site_count
            origin_parent_sites += read_record["cnv_mom_sites"]
            origin_parent_reads += read_record["mom_reads"]
            other_parent_sites += read_record["dad_sites"]
            other_parent_reads += read_record["dad_reads"]
            if ("AMBIGUOUS_READBACKED" in evidence_types):
                evidence_types.remove("AMBIGUOUS_READBACKED")
            evidence_types.append("ALLELE-BALANCE")
    elif ((dad_cnv_site_count + mom_cnv_site_count) > 0) and not ("READBACKED" in evidence_types):
        # this just became ambiguous because of contradictory results
        origin_parent = None
        evidence_count += dad_cnv_site_count + mom_cnv_site_count
        origin_parent_sites += read_record["cnv_dad_sites"]
        other_parent_sites = read_record["cnv_mom_sites"]
        evidence_types.append("AMBIGUOUS_ALLELE-BALANCE")
        ambig = True

    if (origin_parent is None or ambig) and not include_ambiguous:
        return
    origin_parent_sites = sorted(origin_parent_sites)
    other_parent_sites = sorted(other_parent_sites)

    origin_parent_sites = (
        ",".join(origin_parent_sites) if len(origin_parent_sites) > 0 else "-"
    )
    origin_parent_reads = (
        ",".join(origin_parent_reads) if len(origin_parent_reads) > 0 else "-"
    )
    other_parent_sites = (
        ",".join(other_parent_sites) if len(other_parent_sites) > 0 else "-"
    )
    other_parent_reads = (
        ",".join(other_parent_reads) if len(other_parent_reads) > 0 else "-"
    )

    merged_record = {
        "chrom": read_record["region"]["chrom"],
        "start": int(read_record["region"]["start"]),
        "end": int(read_record["region"]["end"]),
        "vartype": read_record["vartype"],
        "kid": read_record["kid"],
        "origin_parent": origin_parent,
        "other_parent": other_parent,
        "evidence_count": evidence_count,
        "evidence_types": evidence_types,
    }
    if verbose:
        merged_record["origin_parent_sites"] = origin_parent_sites
        merged_record["origin_parent_reads"] = origin_parent_reads
        merged_record["other_parent_sites"] = other_parent_sites
        merged_record["other_parent_reads"] = other_parent_reads
    return merged_record


def write_vcf_output(in_vcf_name, read_records, include_ambiguous, verbose, outfile):
    vcf = VCF(in_vcf_name)
    vcf.add_format_to_header(
        {
            "ID": "UOPS",
            "Description": "Count of pieces of evidence supporting the "
            + "unfazed-identified origin parent or `-1` if missing",
            "Type": "Float",
            "Number": "1",
        }
    )
    vcf.add_format_to_header(
        {
            "ID": "UET",
            "Description": "Unfazed evidence type: "
            + "`0` (readbacked), "
            + "`1` (allele-balance, for CNVs only), "
            + "`2` (both), "
            + "`3` (ambiguous readbacked), "
            + "`4` (ambiguous allele-balance), "
            + "`5` (ambiguous both), "
            + "`6` (auto-phased sex-chromosome variant in male), or "
            + "`-1` (missing)",
            "Type": "Float",
            "Number": "1",
        }
    )
    writer = Writer(outfile, vcf)

    for variant in vcf:
        # keys = []
        unfazed_gts = variant.genotypes
        uops = []
        uet = []
        for i, gt in enumerate(variant.gt_types):
            uops_entry = -1
            uet_entry = -1

            if gt in [HET, HOM_ALT]:
                vartype = variant.INFO.get("SVTYPE")
                if vartype is None:
                    vartype = SNV_TYPE

                key_fields = {
                    "chrom": variant.CHROM,
                    "start": variant.start,
                    "end": variant.end,
                    "sample": vcf.samples[i],
                    "vartype": vartype,
                }
                key = "{chrom}_{start}_{end}_{sample}_{vartype}".format(**key_fields)
                if key in read_records:
                    record_summary = summarize_record(
                        read_records[key], include_ambiguous, verbose
                    )
                    if record_summary is not None:
                        origin_parent = record_summary["origin_parent"]
                        if origin_parent == read_records[key]["dad"]:
                            unfazed_gts[i][0] = 1
                            unfazed_gts[i][1] = 0
                            unfazed_gts[i][2] = True
                        elif origin_parent == read_records[key]["mom"]:
                            unfazed_gts[i][0] = 0
                            unfazed_gts[i][1] = 1
                            unfazed_gts[i][2] = True

                        uops_entry = record_summary["evidence_count"]

                        evidence_types = record_summary["evidence_types"]
                        uet_entry = -1
                        if "AMBIGUOUS_READBACKED" in evidence_types:
                            uet_entry = 3
                        elif "AMBIGUOUS_ALLELE-BALANCE" in evidence_types:
                            uet_entry = 4
                        elif "AMBIGUOUS_BOTH" in evidence_types:
                            uet_entry = 5
                        elif "SEX-CHROM" in evidence_types:
                            uet_entry = 6
                        elif (
                            "READBACKED" in evidence_types
                            and "ALLELE-BALANCE" in evidence_types
                        ):
                            uet_entry = 2
                        elif "READBACKED" in evidence_types:
                            uet_entry = 0
                        elif "ALLELE-BALANCE" in evidence_types:
                            uet_entry = 1

            uops.append(uops_entry)
            uet.append(uet_entry)
        variant.genotypes = unfazed_gts
        variant.set_format("UOPS", np.array(uops))
        variant.set_format("UET", np.array(uet))

        writer.write_record(variant)


def write_bed_output(read_records, include_ambiguous, verbose, outfile):
    header = [
        "#chrom",
        "start",
        "end",
        "vartype",
        "kid",
        "origin_parent",
        "other_parent",
        "evidence_count",
        "evidence_types",
    ]

    template_fields = [
        "{chrom}",
        "{start}",
        "{end}",
        "{vartype}",
        "{kid}",
        "{origin_parent}",
        "{other_parent}",
        "{evidence_count}",
        "{evidence_types}",
    ]

    if verbose:
        header += [
            "origin_parent_sites",
            "origin_parent_reads",
            "other_parent_sites",
            "other_parent_reads",
        ]

        template_fields += [
            "{origin_parent_sites}",
            "{origin_parent_reads}",
            "{other_parent_sites}",
            "{other_parent_reads}",
        ]

    template = "\t".join(template_fields)
    record_summaries = []

    for key in read_records:
        record_summary = summarize_record(read_records[key], include_ambiguous, verbose)

        if record_summary is not None:
            record_summaries.append(record_summary)

    record_summaries = sorted(
        record_summaries, key=lambda x: (x["chrom"], x["start"], x["end"])
    )
    if outfile == "/dev/stdout":
        print("\t".join(header))
        for record_summary in record_summaries:
            record_summary["evidence_types"] = ",".join(
                record_summary["evidence_types"]
            )
            print(template.format(**record_summary))
    else:
        with open(outfile, "w") as outfile_fh:
            print("\t".join(header), file=outfile_fh)
            for record_summary in record_summaries:
                record_summary["evidence_types"] = ",".join(
                    record_summary["evidence_types"]
                )
                print(template.format(**record_summary), file=outfile_fh)


def unfazed(args):
    input_type = ""
    bam_names_dict = get_bam_names(args.bam_dir, args.bam_pairs, args.reference)
    snvs = []
    svs = []
    reader = ""
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
        sys.exit(
            "dnms file type is unrecognized. Must be bed, bed.gz, vcf, vcf.gz, or bcf"
        )

    output_type = args.output_type if args.output_type is not None else input_type
    if output_type == "vcf" and input_type != "vcf":
        print(
            "Invalid option: --output-type is vcf, but input is not a vcf type. "
            + "Rerun with `--output-type bed` or input dnms as one of the following:",
            ", ".join(VCF_TYPES),
            file=sys.stderr,
        )
        sys.exit(1)

    kids = set()
    missing_samples = set()
    duplicated_samples = set()
    for var_fields in reader(args.dnms):
        sample = var_fields["kid"]
        if sample not in bam_names_dict:
            if not (sample in missing_samples):
                print("missing alignment file for", sample, file=sys.stderr)
                missing_samples.add(sample)
            continue
        elif len(bam_names_dict[sample]) != 1:
            if not (sample in duplicated_samples):
                print(
                    "multiple alignment files for",
                    sample + ".",
                    "Please specify correct alignment file using --bam-pairs",
                    file=sys.stderr,
                )
                duplicated_samples.add(sample)
            continue
        kids.add(sample)
        bam = list(bam_names_dict[sample])[0]
        var_fields["bam"] = bam
        var_fields["cram_ref"] = args.reference

        if var_fields["vartype"] in SV_TYPES:
            svs.append(var_fields)
        elif var_fields["vartype"] == SNV_TYPE:
            snvs.append(var_fields)

    pedigrees = parse_ped(args.ped, kids)
    kids = list(pedigrees.keys())

    filtered_snvs = []
    for snv in snvs:
        if snv["kid"] in kids:
            filtered_snvs.append(snv)
    snvs = filtered_snvs
    filtered_svs = []
    for sv in svs:
        if sv["kid"] in kids:
            filtered_svs.append(sv)
    svs = filtered_svs

    phased_svs = {}
    phased_snvs = {}

    if (len(snvs) + len(svs)) == 0:
        sys.exit("No phaseable variants")
    if len(svs) > 0:
        phased_svs = phase_svs(
            svs, kids, pedigrees, args.sites, args.threads, args.build, args.no_extended
        )
    if len(snvs) > 0:
        phased_snvs = phase_snvs(
            snvs, kids, pedigrees, args.sites, args.threads, args.build, args.no_extended
        )

    all_phased = phased_snvs
    all_phased.update(phased_svs)

    if output_type == "vcf":
        write_vcf_output(
            args.dnms, all_phased, args.include_ambiguous, args.verbose, args.outfile,
        )
    elif output_type == "bed":
        write_bed_output(all_phased, args.include_ambiguous, args.verbose, args.outfile)
