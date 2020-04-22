#! /usr/bin/env python
from __future__ import print_function
from cyvcf2 import VCF
from concurrent.futures import ThreadPoolExecutor,wait
import sys

HOM_ALT= 3
HET = 1
HOM_REF = 0

SEX_KEY = {
    'male' : 1,
    'female': 2
}

def get_prefix(vcf):
    chrom_prefix = ""
    for var in vcf:
        if "chr" in var.CHROM.lower():
            chrom_prefix = var.CHROM[:3]
        return chrom_prefix
    return ""

def get_position(vcf, denovo, extra, whole_region):
    locs = []
    loc_template = "{prefix}{chrom}:{start}-{end}"
    prefix = get_prefix(vcf)
    if whole_region:
        locs.append(loc_template.format(
            prefix=prefix,
            chrom=denovo['chrom'].strip("chr"), 
            start=(int(denovo['start']) - extra), 
            end=(int(denovo['end']) + extra)
        ))
    else:
        locs.append(loc_template.format(
            prefix=prefix,
            chrom=denovo['chrom'].strip("chr"), 
            start=(int(denovo['start']) - extra), 
            end=(int(denovo['start']) + extra)
        ))
        if (int(denovo['end'])-int(denovo['start'])) > extra:
            locs.append(loc_template.format(
                prefix=prefix,
                chrom=denovo['chrom'].strip("chr"), 
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

    if min_ab <= allele_bal <= max_ab: 
        return True
    return False


def get_kid_allele(denovo, genotypes, ref_depths, alt_depths, kid_idx):
    kid_allele = None
    if ((denovo['svtype'] == "DEL") and (ref_depths[kid_idx] + alt_depths[kid_idx]) > 4):
        #large deletions can be genotyped by hemizygous inheritance of informative alleles
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
        #large duplications can be genotyped by unbalanced het inheritance of informative alleles
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



def find(dnms, pedigrees, vcf_name, search_dist, threads, whole_region=True):
    """
    Given list of denovo variant positions 
    a vcf_name, and the distance upstream or downstream to search, find informative sites
    """
    if len(dnms) > 10000:
        return find_many(dnms, pedigrees, vcf_name, search_dist, threads, whole_region)
    elif len(dnms) <= 0:
        return

    vcf = VCF(vcf_name)
    sample_dict = dict(zip(vcf.samples, range(len(vcf.samples))))
    for i,denovo in enumerate(dnms):
        kid_id = denovo['kid']
        dad_id = pedigrees[denovo['kid']]['dad']
        mom_id = pedigrees[denovo['kid']]['mom']
        missing = False
        for sample_id in [kid_id,dad_id,mom_id]:
            if sample_id not in sample_dict:
                print("{} missing from SNV bcf")
                missing = True
        if missing:
            continue

        try:
            kid_idx = sample_dict[denovo['kid']]
            dad_idx = sample_dict[pedigrees[denovo['kid']]['dad']]
            mom_idx = sample_dict[pedigrees[denovo['kid']]['mom']]
        except Exception as e:
            continue
        
        candidate_sites = [] 
        het_sites = []
        # loop over all variants in the VCF within search_dist bases from the DNM
        for variant in get_position(vcf, denovo, search_dist, whole_region):
            # ignore more complex variants for now
            if (len(variant.ALT) != 1 or (len(variant.REF) > 1) or ('*' in variant.ALT or len(variant.ALT[0]) > 1)): 
                continue

            #male chrX variants have to come from mom
            if variant.CHROM == 'X' and (pedigrees[denovo['kid']]['sex'] == SEX_KEY['male']): 
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
            if ((genotypes[kid_idx] == HET) and 
                    is_high_quality_site(dad_idx, ref_depths, alt_depths, genotypes, gt_quals) and
                    is_high_quality_site(mom_idx, ref_depths, alt_depths, genotypes, gt_quals)):
                #variant usable for extended read-backed phasing
                het_sites.append( {
                    'pos'           : variant.start,
                    'ref_allele'    : variant.REF,
                    'alt_allele'    : variant.ALT[0],
                })


            if whole_region and ('svtype' in denovo):
                candidate['kid_allele'] = get_kid_allele(denovo, genotypes, ref_depths, alt_depths, kid_idx)
                if not candidate['kid_allele']:
                    continue
            elif genotypes[kid_idx] != HET or not is_high_quality_site(kid_idx, ref_depths, alt_depths, genotypes, gt_quals):
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
    
            #if kid is hemizygous we need to make sure the inherited allele is not shared
            #by both parents
            if genotypes[kid_idx] in [HOM_ALT, HOM_REF]:
                unique_allele = True
                #if one parent is het and the other is homozygous for either allele
                #make sure the kid doesn't have that allele
                parent_gts = [genotypes[dad_idx], genotypes[mom_idx]]
                if (HET in parent_gts) and (HOM_ALT in parent_gts or HOM_REF in parent_gts):
                    for parent_gt in parent_gts:
                        kid_gt = genotypes[kid_idx]
                        if (parent_gt in [HOM_ALT,HOM_REF]) and (kid_gt == parent_gt):
                            unique_allele=False
                if not unique_allele:
                    continue
            candidate_sites.append(candidate)
            #if this is an informative site, then don't use it as a het
            #if genotypes[kid_idx] == HET:
            #    het_sites.pop()


        denovo['candidate_sites'] = sorted(candidate_sites, key=lambda x: x['pos'])
        denovo['het_sites'] = sorted(het_sites, key=lambda x: x['pos'])
        dnms[i] = denovo

    return dnms

def create_lookups(dnms):
    #this will be a lookup to find samples for a range where variants are informative for a given denovo
    samples_by_location = {}
    vars_by_sample = {}
    chrom_ranges = {}
    for denovo in dnms:
        chrom = denovo['chrom']
        start = int(denovo['start'])
        end = int(denovo['end'])
        sample = denovo['kid']

        #find the range within the chromosome where we need data
        if chrom not in chrom_ranges:
            chrom_ranges[chrom] = [start,end]
        if start<chrom_ranges[chrom][0]:
            chrom_ranges[chrom][0] = start
        if end>chrom_ranges[chrom][1]:
            chrom_ranges[chrom][1] = end

        if sample not in vars_by_sample:
            vars_by_sample[sample] = {}
        if chrom not in vars_by_sample[sample]:
            vars_by_sample[sample][chrom] = {}
        if start not in vars_by_sample[sample][chrom]:
            vars_by_sample[sample][chrom][start] = []
        vars_by_sample[sample][chrom][start].append(denovo)

        if chrom not in samples_by_location:
            samples_by_location[chrom] = {}

        if start not in samples_by_location[chrom]:
            samples_by_location[chrom][start] = []
        samples_by_location[chrom][start].append(sample)

        if (end-start) > 2:
            if end not in samples_by_location[chrom] and ((end-start) > 2):
                samples_by_location[chrom][end] = []
            samples_by_location[chrom][end].append(sample)
    return samples_by_location,vars_by_sample,chrom_ranges

def get_close_vars(chrom, pos, samples_by_location, vars_by_sample, search_dist, whole_region):
    close_var_keys = []
    if not chrom in samples_by_location:
        return close_var_keys
    if not whole_region:
        for dn_loc in samples_by_location[chrom]:
            if (dn_loc-search_dist) <= pos <= (dn_loc+search_dist):
                close_samples = samples_by_location[chrom][dn_loc]
                for close_sample in close_samples:
                    close_var_keys.append([close_sample,chrom,dn_loc])
    else:
        for dn_loc in samples_by_location[chrom]:
            close_samples = samples_by_location[chrom][dn_loc]
            for close_sample in close_samples:
                for denovo in vars_by_sample[close_sample][chrom][dn_loc]:
                    dn_start = int(denovo["start"])
                    dn_end = int(denovo["end"])
                    if (dn_start-search_dist) <= pos <= (dn_end+search_dist):
                        close_var_keys.append([close_sample,chrom,dn_start])
    return close_var_keys

def get_family_indexes(kid, pedigrees,sample_dict):
    dad_id = pedigrees[kid]['dad']
    mom_id = pedigrees[kid]['mom']
    missing = False
    for sample_id in [kid,dad_id,mom_id]:
        if sample_id not in sample_dict:
            print("{} missing from SNV bcf")
            missing = True
    if missing:
        return None,None,None

    try:
        kid_idx = sample_dict[kid]
        dad_idx = sample_dict[pedigrees[kid]['dad']]
        mom_idx = sample_dict[pedigrees[kid]['mom']]
        return kid_idx,dad_idx,mom_idx
    except Exception as e:
        return None,None,None


def add_good_candidate_variant(variant, vars_by_sample, dn_key, pedigrees, whole_region, sample_dict):
    kid,chrom,pos = dn_key
    kid_idx,dad_idx,mom_idx = get_family_indexes(kid, pedigrees,sample_dict)
    dad = pedigrees[kid]['dad']
    mom = pedigrees[kid]['mom']
    if None in [kid_idx,dad_idx,mom_idx]:
        return False

    for i,denovo in enumerate(vars_by_sample[kid][chrom][pos]):
        #male chrX variants have to come from mom
        if variant.CHROM == 'X' and (pedigrees[kid]['sex'] == SEX_KEY['male']): 
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
        if ((genotypes[kid_idx] == HET) and 
                is_high_quality_site(dad_idx, ref_depths, alt_depths, genotypes, gt_quals) and
                is_high_quality_site(mom_idx, ref_depths, alt_depths, genotypes, gt_quals)):
            if not "het_sites" in vars_by_sample[kid][chrom][pos][i]:
                vars_by_sample[kid][chrom][pos][i]["het_sites"] = []
            #variant usable for extended read-backed phasing
            vars_by_sample[kid][chrom][pos][i]["het_sites"].append( {
                'pos'           : variant.start,
                'ref_allele'    : variant.REF,
                'alt_allele'    : variant.ALT[0],
            })


        if whole_region and ('svtype' in denovo):
            candidate['kid_allele'] = get_kid_allele(denovo, genotypes, ref_depths, alt_depths, kid_idx)
            if not candidate['kid_allele']:
                continue
        elif genotypes[kid_idx] != HET or not is_high_quality_site(kid_idx, ref_depths, alt_depths, genotypes, gt_quals):
            continue

        if not (is_high_quality_site(dad_idx, ref_depths, alt_depths, genotypes, gt_quals) and  
                is_high_quality_site(mom_idx, ref_depths, alt_depths, genotypes, gt_quals)): 
            continue

        if genotypes[dad_idx] in (HET, HOM_ALT) and genotypes[mom_idx] == HOM_REF:
            candidate['alt_parent'] = dad
            candidate['ref_parent'] = mom
        elif genotypes[mom_idx] in (HET, HOM_ALT) and genotypes[dad_idx] == HOM_REF:
            candidate['alt_parent'] = mom
            candidate['ref_parent'] = dad
        elif genotypes[mom_idx] == HET and genotypes[dad_idx] == HOM_ALT:
            candidate['alt_parent'] = dad
            candidate['ref_parent'] = mom
        elif genotypes[dad_idx] == HET and genotypes[mom_idx] == HOM_ALT:
            candidate['alt_parent'] = mom
            candidate['ref_parent'] = dad
        else:
            continue

        #if kid is hemizygous we need to make sure the inherited allele is not shared
        #by both parents
        if genotypes[kid_idx] in [HOM_ALT, HOM_REF]:
            unique_allele = True
            #if one parent is het and the other is homozygous for either allele
            #make sure the kid doesn't have that allele
            parent_gts = [genotypes[dad_idx], genotypes[mom_idx]]
            if (HET in parent_gts) and (HOM_ALT in parent_gts or HOM_REF in parent_gts):
                for parent_gt in parent_gts:
                    kid_gt = genotypes[kid_idx]
                    if (parent_gt in [HOM_ALT,HOM_REF]) and (kid_gt == parent_gt):
                        unique_allele=False
            if not unique_allele:
                continue
        if not "candidate_sites" in vars_by_sample[kid][chrom][pos][i]:
            vars_by_sample[kid][chrom][pos][i]['candidate_sites'] = []
        vars_by_sample[kid][chrom][pos][i]['candidate_sites'].append(candidate)
    return True

#########################################################################################################################
def multithread_find_many(vcf_name,chrom, chrom_range,samples_by_location, vars_by_sample, search_dist, pedigrees, whole_region):
    vcf = VCF(vcf_name)
    prefix = get_prefix(vcf)
    sample_dict = dict(zip(vcf.samples, range(len(vcf.samples))))
    search_string = "{}:{}-{}".format(prefix+chrom.strip("chr"), chrom_range[0]-search_dist, chrom_range[1]+search_dist)
    vcf_iter = vcf(search_string)
    
    for variant in vcf_iter:

        # ignore complex variants for now
        if (len(variant.ALT) != 1 or (len(variant.REF) > 1) or ('*' in variant.ALT or len(variant.ALT[0]) > 1)): 
            continue

        close_var_keys = get_close_vars(variant.CHROM, variant.POS, samples_by_location, vars_by_sample, search_dist, whole_region)

        if len(close_var_keys) == 0:
            continue
        for close_var_key in close_var_keys:
            add_good_candidate_variant(variant, vars_by_sample, close_var_key, pedigrees, whole_region, sample_dict)



def find_many(dnms, pedigrees, vcf_name, search_dist, threads, whole_region=True):
    """
    Given list of denovo variant positions 
    a vcf_name, and the distance upstream or downstream to search, find informative sites
    """
    samples_by_location,vars_by_sample,chrom_ranges = create_lookups(dnms)
    chroms = set([dnm['chrom'] for dnm in dnms])

    executor = ThreadPoolExecutor(threads)
    futures = []
    for chrom in chroms:
        futures.append(executor.submit(
            multithread_find_many, 
            vcf_name,
            chrom, 
            chrom_ranges[chrom],
            samples_by_location, 
            vars_by_sample, 
            search_dist, 
            pedigrees,
            whole_region
        ))
    
    wait(futures)

    dnms_annotated = []
    for sample in vars_by_sample:
        for chrom in vars_by_sample[sample]:
            for pos in vars_by_sample[sample][chrom]:
                for denovo in vars_by_sample[sample][chrom][pos]:
                    if 'candidate_sites' in denovo:
                        denovo['candidate_sites'] = sorted(denovo['candidate_sites'], key=lambda x: x['pos'])
                    if 'het_sites' in denovo:
                        denovo['het_sites'] = sorted(denovo['het_sites'], key=lambda x: x['pos'])
                    dnms_annotated.append(denovo)

    return dnms_annotated


if __name__ == "__main__":
    sys.exit("Import this as a module")
