#! /usr/bin/env python
from __future__ import print_function
import argparse
import toolshed as ts
from cyvcf2 import VCF

def get_position(vcf, d, ref=None, alt=None, extra=0):
    loc = "%s:%d-%d" % (d['chrom'], int(d['end']) - extra,
                        int(d['end']) + extra)
    for v in vcf(loc):
        if ref and v.REF != ref: continue
        if alt and alt not in v.ALT: continue
        yield v

def is_high_quality_site(i, rd, ad, gts, gqs, 
                            min_gq=20, min_depth=10):
    """
    check if the potential informative site is high quality.
    this is pretty hacky, but filters out the really cruddy variants.
    i: index of the parent in the VCF FORMAT field
    rd: numpy array of reference depths
    ad: numpy array of alternate depths
    gts: numpy array of genotypes
    gqs: numpy array of genotype qualities
    """
    HOM_ALT, HET = 3, 1
    if gts[i] == HOM_ALT:
        min_ab, max_ab = 0.8, 1.
    elif gts[i] == HET:
        min_ab, max_ab = 0.2, 0.8
    if gqs[i] < min_gq: 
        return False
    if rd[i] + ad[i] < min_depth:
        return False
    if ad[i] / float(rd[i] + ad[i]) < min_ab: 
        return False
    if ad[i] / float(rd[i] + ad[i]) > max_ab: 
        return False
    return True

def run(args, fragment_length=500):

    HOM_ALT, HET, HOM_REF = 3, 1, 0

    vcf = VCF(args.vcf)

    sample_dict = dict(zip(vcf.samples, range(len(vcf.samples))))
    import gzip
    import sys
    import os
    import time
    dnms = gzip.open(args.dnms) if args.dnms.endswith('.gz') else args.dnms
    for i, d in enumerate(ts.reader(args.dnms, header="ordered")):
        if i>0 and i % 10 == 0: 
            print ("done with {} variants".format(i),file=sys.stderr)
        if i == 0: 
            keys = list(d.keys())
            keys.extend(['candidate_sites'])
            print('\t'.join(keys))
        if args.chrom and d['chrom'] != args.chrom: continue
        ikid = sample_dict[d['sample_id']]
        idad = sample_dict[d['paternal_id']]
        imom = sample_dict[d['maternal_id']]
        
        candidate_sites = [] 

        # loop over all variants in the VCF +/- 350 bp from the DNM
        t = time.time()
        for v in get_position(vcf, d, extra=fragment_length):
            if len(candidate_sites) > 3: break 
            # ignore more complex variants for now
            if len(v.REF) > 1: continue
            if len(v.ALT) > 1 or '*' in v.ALT or len(v.ALT[0]) > 1: continue
            gts = v.gt_types
            # special logic for male chrX (hemizygous) variants
            is_male = d['sample_sex'] == 'male'
            if v.CHROM == 'X' and is_male and gts[ikid] != HOM_ALT: continue
            else:
                if gts[ikid] != HET: continue
            candidate_parent = 'NA'
            candidate_site = 'NA'
            rd, ad = v.gt_ref_depths, v.gt_alt_depths
            gqs = v.gt_quals
            if gts[idad] in (HET, HOM_ALT) and gts[imom] == HOM_REF:
                if not is_high_quality_site(idad, rd, ad, gts, gqs): continue
                candidate_site = str(v.start) + ":" + str(v.REF) + '-' + str(v.ALT[0])
                candidate_parent = d['paternal_id']

            elif gts[imom] in (HET, HOM_ALT) and gts[idad] == HOM_REF:
                if not is_high_quality_site(imom, rd, ad, gts, gqs): continue
                candidate_site = str(v.start) + ":" + str(v.REF) + '-' + str(v.ALT[0])
                candidate_parent = d['maternal_id']
    
            candidate_sites.append(candidate_parent + '|' + candidate_site)
        val = ','.join([x for x in candidate_sites if 'NA' not in x])
        if len(val) == 0: val = 'NA'
        d['candidate_sites'] = val

        print('\t'.join(d.values()))

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dnms", help="Path to BED files of DNMs to phase", default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/phaser/test/data/dnms.bed")
    p.add_argument("--vcf", help="Path to VCF we'll use to annotate DNMs", default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/data/ceph.bcf")
    p.add_argument("-chrom")
    args = p.parse_args()
    run(args)

