#! /usr/bin/env python
from __future__ import print_function

import sys
from concurrent.futures import ThreadPoolExecutor, wait

from .informative_site_finder import find
from .read_collector import collect_reads_sv
from .site_searcher import match_informative_sites

SEX_KEY = {"male": 1, "female": 2}
QUIET_MODE = False
CONCORDANT_UPPER_LENS = {}
# https://www.ncbi.nlm.nih.gov/grc/human
grch37_par1 = {
    "x": [10001, 2781479],
    "y": [10001, 2781479],
}
grch37_par2 = {
    "x": [155701383, 156030895],
    "y": [56887903, 57217415],
}

grch38_par1 = {
    "x": [60001, 2699520],
    "y": [10001, 2649520],
}

grch38_par2 = {
    "x": [154931044, 155260560],
    "y": [59034050, 59363566],
}


def phase_by_reads(matches):
    # parent_ids -> list of informative site matches
    origin_parent_data = {}

    for ref_alt in matches:
        for match_info in matches[ref_alt]:
            read = match_info["read"]
            for match in match_info["matches"]:
                if len(origin_parent_data) == 0:
                    origin_parent_data[match["ref_parent"]] = []
                    origin_parent_data[match["alt_parent"]] = []
                # to avoid issues from indels, use the reference position to index the read
                try:
                    read_pos = read.get_reference_positions(full_length=True).index(
                        match["pos"]
                    )
                except ValueError:
                    continue
                kid_allele = read.query_sequence[read_pos]

                # if the kid base is ref, this read comes from ref_parent
                # if the kid base is alt, this read comes from alt_parent
                # reads are already assigned to ref or alt relative to the de novo
                # kid base is ref or alt relative to the informative site
                read_origin = ""
                if kid_allele == match["ref_allele"]:
                    read_origin = "ref_parent"
                elif kid_allele == match["alt_allele"]:
                    read_origin = "alt_parent"
                else:
                    continue

                # ref_parent means the parent that has the reference allele at the informative site
                # ref haplotype means the read comes from the non-denovo haplotye
                # so if the read comes from the ref_parent and the ref haplotype,
                # de novo is on the alt parent
                if read_origin == "ref_parent":
                    if ref_alt == "ref":
                        origin_parent_data[match["alt_parent"]].append(
                            [read, match["pos"]]
                        )
                    else:
                        origin_parent_data[match["ref_parent"]].append(
                            [read, match["pos"]]
                        )
                else:
                    if ref_alt == "ref":
                        origin_parent_data[match["ref_parent"]].append(
                            [read, match["pos"]]
                        )
                    else:
                        origin_parent_data[match["alt_parent"]].append(
                            [read, match["pos"]]
                        )
    return origin_parent_data


def phase_by_snvs(informative_sites):
    if len(informative_sites) <= 0:
        return None
    # parent_ids -> list of informative site matches
    origin_parent_data = {
        informative_sites[0]["ref_parent"]: [],
        informative_sites[0]["alt_parent"]: [],
    }

    # split informative sites up by the parent they say is responsible
    for informative_site in informative_sites:
        origin_parent_data[informative_site[informative_site["kid_allele"]]].append(
            informative_site
        )
    return origin_parent_data


def multithread_read_phasing(
    denovo,
    records,
    dad_id,
    mom_id,
    no_extended,
    insert_size_max_sample,
    stdevs,
    min_map_qual,
    min_gt_qual,
    readlen,
    split_error_margin,
):
    region = {
        "chrom": denovo["chrom"],
        "start": denovo["start"],
        "end": denovo["end"],
    }
    concordant_upper_len = None
    if denovo["kid"] in CONCORDANT_UPPER_LENS:
        concordant_upper_len = CONCORDANT_UPPER_LENS[denovo["kid"]]

    # these are reads that support the ref or alt allele of the de novo variant
    dnm_reads,concordant_upper_len = collect_reads_sv(
        denovo["bam"],
        region,
        denovo["het_sites"],
        denovo["cram_ref"],
        no_extended,
        concordant_upper_len,
        insert_size_max_sample,
        stdevs,
        min_map_qual,
        min_gt_qual,
        readlen,
        split_error_margin,
    )
    CONCORDANT_UPPER_LENS[denovo["kid"]] = concordant_upper_len
    matches = match_informative_sites(dnm_reads, denovo["candidate_sites"])

    if len(matches["alt"]) <= 0 and len(matches["ref"]) <= 0:
        if not QUIET_MODE:
            print(
                "No reads overlap informative sites for variant {chrom}:{start}-{end}".format(
                    **region
                ),
                file=sys.stderr,
            )
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
        "region": region,
        "vartype": denovo["vartype"],
        "kid": denovo["kid"],
        "dad": dad_id,
        "mom": mom_id,
        "dad_sites": dad_informative_sites,
        "mom_sites": mom_informative_sites,
        "evidence_type": "readbacked",
        "dad_reads": dad_reads,
        "mom_reads": mom_reads,
        "cnv_dad_sites": "",
        "cnv_mom_sites": "",
        "cnv_evidence_type": "",
    }
    key = [str(v) for v in region.values()] + [denovo["kid"], denovo["vartype"]]
    records["_".join(key)] = record


def run_read_phasing(
    dnms,
    pedigrees,
    vcf,
    threads,
    build,
    no_extended,
    multiread_proc_min,
    ab_homref,
    ab_homalt,
    ab_het,
    min_gt_qual,
    min_depth,
    search_dist,
    insert_size_max_sample,
    stdevs,
    min_map_qual,
    readlen,
    split_error_margin,
):
    # get informative sites near the breakpoints of SVs for reab-backed phasing
    dnms_with_informative_sites = find(
        dnms,
        pedigrees,
        vcf,
        search_dist,
        threads,
        build,
        multiread_proc_min,
        QUIET_MODE,
        ab_homref,
        ab_homalt,
        ab_het,
        min_gt_qual,
        min_depth,
        whole_region=False,
    )
    records = {}
    if threads != 1:
        executor = ThreadPoolExecutor(threads)
        futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo["kid"]]["dad"]
        mom_id = pedigrees[denovo["kid"]]["mom"]
        if autophase(denovo, pedigrees, records, dad_id, mom_id, build):
            continue

        if "candidate_sites" not in denovo or len(denovo["candidate_sites"]) == 0:
            if not QUIET_MODE:
                print(
                    "No usable informative sites for read-based phasing of variant {}:{}-{}".format(
                        denovo["chrom"], denovo["start"], denovo["end"]
                    ),
                    file=sys.stderr,
                )
            continue
        if threads != 1:
            futures.append(
                executor.submit(
                    multithread_read_phasing,
                    denovo,
                    records,
                    dad_id,
                    mom_id,
                    no_extended,
                    insert_size_max_sample,
                    stdevs,
                    min_map_qual,
                    min_gt_qual,
                    readlen,
                    split_error_margin,
                )
            )
        else:
            multithread_read_phasing(
                denovo,
                records,
                dad_id,
                mom_id,
                no_extended,
                insert_size_max_sample,
                stdevs,
                min_map_qual,
                min_gt_qual,
                readlen,
                split_error_margin,
            )
    if threads != 1:
        wait(futures)
    return records


def multithread_cnv_phasing(denovo, records, dad_id, mom_id):
    region = {
        "chrom": denovo["chrom"],
        "start": denovo["start"],
        "end": denovo["end"],
    }

    origin_data = phase_by_snvs(denovo["candidate_sites"])
    if not origin_data:
        return

    evidence_items = {dad_id: [], mom_id: []}
    for parent in evidence_items:
        if parent in origin_data and len(origin_data[parent]) > 0:
            evidence_items[parent] = [str(o["pos"]) for o in origin_data[parent]]

    record = {
        "region": region,
        "vartype": denovo["vartype"],
        "kid": denovo["kid"],
        "dad": dad_id,
        "mom": mom_id,
        "cnv_dad_sites": evidence_items[dad_id],
        "cnv_mom_sites": evidence_items[mom_id],
        "cnv_evidence_type": "ALLELE-BALANCE",
        "dad_sites": "",
        "mom_sites": "",
        "evidence_type": "",
        "dad_reads": [],
        "mom_reads": [],
    }
    key = [str(r) for r in region.values()] + [denovo["kid"], denovo["vartype"]]
    records["_".join(key)] = record


def autophase(denovo, pedigrees, records, dad_id, mom_id, build):
    """
    variants in males on the X or Y chromosome and not in the
    pseudoautosomal regions can be automatically phased to the
    dad (if Y) or mom (if x)
    """
    chrom = denovo["chrom"].lower().strip("chr")
    if chrom not in ["y", "x"]:
        return False
    if int(pedigrees[denovo["kid"]]["sex"]) != SEX_KEY["male"]:
        return False
    if build not in ["37", "38"]:
        return False

    if build == "37":
        par1 = grch37_par1
        par2 = grch37_par2

    if build == "38":
        par1 = grch38_par1
        par2 = grch38_par2
    # variant is pseudoautosomal
    if (
        par1[chrom][0] <= denovo["start"] <= par1[chrom][1]
        or par2[chrom][0] <= denovo["start"] <= par2[chrom][1]
    ):
        return False

    region = {
        "chrom": denovo["chrom"],
        "start": denovo["start"],
        "end": denovo["end"],
    }
    record = {
        "region": region,
        "vartype": denovo["vartype"],
        "kid": denovo["kid"],
        "dad": dad_id,
        "mom": mom_id,
        "cnv_dad_sites": "NA",
        "cnv_mom_sites": "NA",
        "cnv_evidence_type": "SEX-CHROM",
        "dad_sites": "",
        "mom_sites": "",
        "evidence_type": "SEX-CHROM",
        "dad_reads": [],
        "mom_reads": [],
    }

    key = [str(r) for r in region.values()] + [denovo["kid"], denovo["vartype"]]
    records["_".join(key)] = record


def run_cnv_phasing(
    dnms,
    pedigrees,
    vcf,
    threads,
    build,
    multithread_proc_min,
    ab_homref,
    ab_homalt,
    ab_het,
    min_gt_qual,
    min_depth,
):
    """
    Specialized phasing for CNVs,
    using the informative sites from the region with a copy-number change
    """
    # get informative sites inside CNVs for purely SNV-based phasing
    dnms_with_informative_sites = find(
        dnms,
        pedigrees,
        vcf,
        0,
        threads,
        build,
        multithread_proc_min,
        QUIET_MODE,
        ab_homref,
        ab_homalt,
        ab_het,
        min_gt_qual,
        min_depth,
    )
    records = {}
    if threads != 1:
        executor = ThreadPoolExecutor(threads)
        futures = []

    for denovo in dnms_with_informative_sites:
        dad_id = pedigrees[denovo["kid"]]["dad"]
        mom_id = pedigrees[denovo["kid"]]["mom"]
        if autophase(denovo, pedigrees, records, dad_id, mom_id, build):
            continue

        if denovo["vartype"] not in ["DEL", "DUP"]:
            continue

        if "candidate_sites" not in denovo or len(denovo["candidate_sites"]) == 0:
            if not QUIET_MODE:
                print(
                    "No usable informative sites for allele-balance phasing of variant {}:{}-{}".format(
                        denovo["chrom"], denovo["start"], denovo["end"]
                    ),
                    file=sys.stderr,
                )
            continue
        if threads != 1:
            futures.append(
                executor.submit(
                    multithread_cnv_phasing, denovo, records, dad_id, mom_id
                )
            )
        else:
            multithread_cnv_phasing(denovo, records, dad_id, mom_id)
    if threads != 1:
        wait(futures)
    return records


# def phase_svs(args):
def phase_svs(
    dnms,
    kids,
    pedigrees,
    sites,
    threads,
    build,
    no_extended,
    multiread_proc_min,
    quiet_mode,
    ab_homref,
    ab_homalt,
    ab_het,
    min_gt_qual,
    min_depth,
    search_dist,
    insert_size_max_sample,
    stdevs,
    min_map_qual,
    readlen,
    split_error_margin,
):
    global QUIET_MODE
    QUIET_MODE = quiet_mode
    cnv_records = run_cnv_phasing(
        dnms,
        pedigrees,
        sites,
        threads,
        build,
        multiread_proc_min,
        ab_homref,
        ab_homalt,
        ab_het,
        min_gt_qual,
        min_depth,
    )
    read_records = run_read_phasing(
        dnms,
        pedigrees,
        sites,
        threads,
        build,
        no_extended,
        multiread_proc_min,
        ab_homref,
        ab_homalt,
        ab_het,
        min_gt_qual,
        min_depth,
        search_dist,
        insert_size_max_sample,
        stdevs,
        min_map_qual,
        readlen,
        split_error_margin,
    )
    for key in cnv_records:
        if key not in read_records:
            read_records[key] = cnv_records[key]
        else:
            read_records[key]["cnv_dad_sites"] = cnv_records[key]["cnv_dad_sites"]
            read_records[key]["cnv_mom_sites"] = cnv_records[key]["cnv_mom_sites"]
            read_records[key]["evidence_type"] += (
                "," + cnv_records[key]["cnv_evidence_type"]
            )
    return read_records
