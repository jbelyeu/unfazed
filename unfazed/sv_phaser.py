#! /usr/bin/env python
from __future__ import print_function

import sys
from .informative_site_finder import find

def run_cnv_phasing(dnms, pedigrees, vcf, threads, build, multithread_proc_min):
    """
    Specialized phasing for CNVs,
    using the informative sites from the region with a copy-number change
    """
    # get informative sites inside CNVs for purely SNV-based phasing
    dnms_with_informative_sites = find(dnms, pedigrees, vcf, 0, threads, build, multithread_proc_min)
    for dnm in dnms_with_informative_sites:
        if 'candidate_sites' in dnm and len(dnm['candidate_sites']) == 0:
            continue
        print("{}:{}-{}".format(dnm["chrom"], dnm["start"], dnm["end"]))
        sites = []

        for candidate in dnm['candidate_sites']:
            sites.append(candidate['pos'])
            print("{}\t{}\t{}".format(dnm["chrom"],candidate["pos"], candidate["pos"]+1))
        print()

# def phase_svs(args):
def phase_svs(dnms, kids, pedigrees, sites, threads, build, no_extended, multiread_proc_min):
    run_cnv_phasing(dnms, pedigrees, sites, threads, build, multiread_proc_min)
