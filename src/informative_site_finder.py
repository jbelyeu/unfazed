#! /usr/bin/env python
from __future__ import print_function
from cyvcf2 import VCF
import sys


def get_position(vcf, denovo, extra, whole_region):
    locs = []
    loc_template = "{chrom}:{start}-{end}"
    if whole_region:
        locs.append(loc_template.format(
            chrom=denovo['chrom'], 
            start=(int(denovo['start']) - extra), 
            end=(int(denovo['end']) + extra)
        ))
    else:
        locs.append(loc_template.format(
            chrom=denovo['chrom'], 
            start=(int(denovo['start']) - extra), 
            end=(int(denovo['start']) + extra)
        ))
        locs.append(loc_template.format(
            chrom=denovo['chrom'], 
            start=(int(denovo['end']) - extra), 
            end=(int(denovo['end']) + extra)
        ))
    for loc in locs:
        for variant in vcf(loc):
            yield variant
    
def is_high_quality_site(i, ref_depths, alt_depths, genotypes, gt_quals, 
                            min_gq=20, min_depth=10):
    """
    check if the potential informative site is high quality.
    this is pretty hacky, but filters out the really cruddy variants.
    i: index of the parent in the VCF FORMAT field
    ref_depths: numpy array of reference depths
    alt_depths: numpy array of alternate depths
    genotypes: numpy array of genotypes
    gt_quals: numpy array of genotype qualities
    """
    HOM_ALT, HET = 3, 1
    if genotypes[i] == HOM_ALT:
        min_ab, max_ab = 0.8, 1.
    elif genotypes[i] == HET:
        min_ab, max_ab = 0.2, 0.8
    if gt_quals[i] < min_gq: 
        return False
    if ref_depths[i] + alt_depths[i] < min_depth:
        return False
    if alt_depths[i] / float(ref_depths[i] + alt_depths[i]) < min_ab: 
        return False
    if alt_depths[i] / float(ref_depths[i] + alt_depths[i]) > max_ab: 
        return False
    return True

def find(dnms, pedigrees, vcf_name, search_dist, whole_region=True):
    """
    Given list of denovo variant positions 
    a vcf_name, and the distance upstream or downstream to search, find informative sites
    """
    
    if len(dnms) <= 0:
        return
    
    HOM_ALT, HET, HOM_REF = 3, 1, 0
    SEX_KEY = {
        'male' : 1,
        'female': 2
    }
    

    vcf = VCF(vcf_name)
    sample_dict = dict(zip(vcf.samples, range(len(vcf.samples))))
    for i,denovo in enumerate(dnms):
        kid_idx = sample_dict[denovo['kid']]
        dad_idx = sample_dict[pedigrees[denovo['kid']]['dad']]
        mom_idx = sample_dict[pedigrees[denovo['kid']]['mom']]
        
        candidate_sites = [] 

        # loop over all variants in the VCF within search_dist bases from the DNM
        for variant in get_position(vcf, denovo, search_dist, whole_region):
            # ignore more complex variants for now
            if ((len(variant.REF) > 1) or 
                        (len(variant.ALT) > 1) or 
                        ('*' in variant.ALT or len(variant.ALT[0]) > 1)): 
                continue

            #male chrX variants have to come from mom
            if variant.CHROM == 'X' and (denovo['sex'] == SEX_KEY['male']): 
                continue

            genotypes = variant.gt_types
            #if the kid's genotype isn't het it can't be used as an informative site
            if genotypes[kid_idx] != HET: continue

            candidate_parent = 'NA'
            candidate_site = 'NA'
            ref_depths = variant.gt_ref_depths
            alt_depths = variant.gt_alt_depths
            gt_quals = variant.gt_quals

            candidate = {
                'pos'           : variant.start,
                'ref_allele'    : variant.REF,
                'alt_allele'    : variant.ALT[0],
            }
            
            if genotypes[dad_idx] in (HET, HOM_ALT) and genotypes[mom_idx] == HOM_REF:
                if not is_high_quality_site(dad_idx, ref_depths, alt_depths, genotypes, gt_quals): 
                    continue
                candidate['alt_parent'] = pedigrees[denovo['kid']]['dad']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['mom']
                candidate_sites.append(candidate)

            elif genotypes[mom_idx] in (HET, HOM_ALT) and genotypes[dad_idx] == HOM_REF:
                if not is_high_quality_site(mom_idx, ref_depths, alt_depths, genotypes, gt_quals): 
                    continue
                candidate['alt_parent'] = pedigrees[denovo['kid']]['mom']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['dad']
                candidate_sites.append(candidate)
        denovo['candidate_sites'] = sorted(candidate_sites, key=lambda x: x['pos'])
        dnms[i] = denovo
    return dnms

if __name__ == "__main__":
    sys.exit("Import this as a module")
