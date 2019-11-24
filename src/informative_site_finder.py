#! /usr/bin/env python
from __future__ import print_function
from cyvcf2 import VCF
import sys

HOM_ALT= 3
HET = 1
HOM_REF = 0

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
    if genotypes[i] == HOM_REF:
        min_ab, max_ab = 0.0, 0.2
    elif genotypes[i] == HOM_ALT:
        min_ab, max_ab = 0.8, 1.
    elif genotypes[i] == HET:
        min_ab, max_ab = 0.2, 0.8
    if gt_quals[i] < min_gq: 
        return False
    if (ref_depths[i] + alt_depths[i]) < min_depth:
        return False

    allele_bal = float(alt_depths[i] / float(ref_depths[i] + alt_depths[i]))
    if allele_bal < min_ab: 
        return False
    if allele_bal > max_ab: 
        return False
    return True


def get_kid_allele(denovo, genotypes, ref_depths, alt_depths, kid_idx):
    kid_allele = None
    if ((denovo['svtype'] == "DEL") and (ref_depths[kid_idx] + alt_depths[kid_idx]) > 4):
        #large deletions can be gentyped by hemizygous inheritance of informative alleles
        if (genotypes[kid_idx] == HOM_ALT):
            kid_allele = 'ref_parent'   
        elif (genotypes[kid_idx] == HOM_REF):
            kid_allele = 'alt_parent'
        else:
            #het kid, variant unusable
            return
    elif ((denovo['svtype'] == "DUP") 
            and (ref_depths[kid_idx] > 2)
            and (alt_depths[kid_idx] > 2)
            and (ref_depths[kid_idx] + alt_depths[kid_idx]) > 10):
        #large duplications can be gentyped by unbalanced het inheritance of informative alleles
        #if there's enough depth
        if (genotypes[kid_idx] == HET):
            kid_alt_allele_bal = (alt_depths[kid_idx]/float(ref_depths[kid_idx]+alt_depths[kid_idx]))
            #allele balance should be about 2:1 for the dominant allele
            if (kid_alt_allele_bal >= .67):
                #this variant came from the alt parent
                kid_allele = 'alt_parent'
            elif (kid_alt_allele_bal <= .33):
                #this variant came from the ref parent
                kid_allele = 'ref_parent'
            else:
                #balanced allele, variant unusable
                return
        else:
            #homozygous kid, variant unusable
            return
    else:
        #not a CNV, or not enough depth
        return
    return kid_allele


def find(dnms, pedigrees, vcf_name, search_dist, whole_region=True):
    """
    Given list of denovo variant positions 
    a vcf_name, and the distance upstream or downstream to search, find informative sites
    """
    
    if len(dnms) <= 0:
        return
    
    SEX_KEY = {
        'male' : 1,
        'female': 2
    }
    #for variant-phasing rather than read-backed, dels require hemizygous and dups require hets
    SV_INFORMATIVE_KEY = {
        'DUP' : HET,
        'DEL' : HOM_ALT
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
            if ((len(variant.REF) > 1) or (len(variant.ALT) > 1) or ('*' in variant.ALT or len(variant.ALT[0]) > 1)): 
                continue

            #male chrX variants have to come from mom
            if variant.CHROM == 'X' and (denovo['sex'] == SEX_KEY['male']): 
                continue

            genotypes = variant.gt_types
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


            if whole_region and ('svtype' in denovo):
                candidate['kid_allele'] = get_kid_allele(denovo, genotypes, ref_depths, alt_depths, kid_idx)
                if not candidate['kid_allele']:
                    continue
            elif genotypes[kid_idx] != HET:
                continue

            if not (is_high_quality_site(dad_idx, ref_depths, alt_depths, genotypes, gt_quals) and  
                    is_high_quality_site(mom_idx, ref_depths, alt_depths, genotypes, gt_quals)): 
                continue

            if genotypes[dad_idx] in (HET, HOM_ALT) and genotypes[mom_idx] == HOM_REF:
                candidate['alt_parent'] = pedigrees[denovo['kid']]['dad']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['mom']
            elif genotypes[mom_idx] in (HET, HOM_ALT) and genotypes[dad_idx] == HOM_REF:
                candidate['alt_parent'] = pedigrees[denovo['kid']]['mom']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['dad']
            elif genotypes[mom_idx] == HET and genotypes[dad_idx] == HOM_ALT:
                candidate['alt_parent'] = pedigrees[denovo['kid']]['dad']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['mom']
            elif genotypes[dad_idx] == HET and genotypes[mom_idx] == HOM_ALT:
                candidate['alt_parent'] = pedigrees[denovo['kid']]['mom']
                candidate['ref_parent'] = pedigrees[denovo['kid']]['dad']
            else:
                continue
    
            #if kid is hemizygous we need to make sure the allele is not shared
            #by both parents
            if genotypes[kid_idx] in [HOM_ALT, HOM_REF]:
                unique_allele = True
                for parent_idx in dad_idx,mom_idx:
                    parent_gt = genotypes[parent_idx]
                    kid_gt = genotypes[kid_idx]
                    if (parent_gt in [HOM_ALT,HOM_REF]) and (genotypes[kid_idx] == parent_gt):
                        unique_allele=False
                if not unique_allele:
                    continue
            candidate_sites.append(candidate)
        denovo['candidate_sites'] = sorted(candidate_sites, key=lambda x: x['pos'])
        dnms[i] = denovo
    return dnms

if __name__ == "__main__":
    sys.exit("Import this as a module")
