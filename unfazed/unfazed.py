#!/usr/bin/env python
from __future__ import print_function

import gzip
import os
import sys
from glob import glob

import numpy as np
from cyvcf2 import VCF, Writer

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

def unfazed(args):
    input_type = ""
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
        kids.add(sample)
        var_fields["cram_ref"] = args.reference

        if var_fields["vartype"] in SV_TYPES:
            svs.append(var_fields)

    pedigrees = parse_ped(args.ped, kids)
    kids = list(pedigrees.keys())

    filtered_svs = []
    for sv in svs:
        if sv["kid"] in kids:
            filtered_svs.append(sv)
    svs = filtered_svs

    phased_svs = {}

    if (len(svs)) == 0:
        sys.exit("No phaseable variants")
    if len(svs) > 0:
        phased_svs = phase_svs(
            svs, kids, pedigrees, args.sites, args.threads, args.build, args.no_extended, args.multiread_proc_min
        )
