from __future__ import absolute_import, print_function
import time
import sys
import os
import pysam
import toolshed as ts
from peddy import Ped
from collections import defaultdict

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("dnms", help="bed file of denovos with a column header for 'chrom,start,end,sample_id")
    p.add_argument("ped", help="ped file with paths to CRAM files for each sample")
    p.add_argument("-chrom")
    args = p.parse_args(argv)
    run(args)

def indel_in_read(read):
    cts = read.cigartuples

    if any([x[0] in (1, 2) for x in cts]): return True
    else:return False

def clips_at_read_start(read):
    """
    ct: pysam cigartuples() object
    """
    SOFT_CLIP = 4
    HARD_CLIP = 5

    cts = read.cigartuples
    if cts[0][0] in (SOFT_CLIP, HARD_CLIP):
        return cts[0][1]
    else: return 0

def get_dnm_base(cram, chrom, pos, read=''):

    dnm_base = None 
    for pileupcolumn in cram.pileup(chrom, pos - 1, pos + 1):
        # do everything w/r/t the informative site
        if pileupcolumn.pos < pos: continue
        if pileupcolumn.pos > pos: break
        # for every read that spans the informative site, check that 
        # the read has support for both the informative HET and the known DNM.
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.query_name != read: continue
            read_seq = pileupread.alignment.query_sequence
            read_pos = pileupread.query_position
            if pileupread.is_del or pileupread.is_refskip: continue
            if pileupread.indel > 0:
                dnm_base = '+'
            elif pileupread.indel < 0:
                dnm_base = '-'
            else:
                try:
                    dnm_base = read_seq[read_pos]
                except IndexError: continue

    return dnm_base

def phased_to_inf(rd, inf_read_base='A', dnm_read_base='T'):
    """
    check if the DNM is phased with HET variant in the informative parent,
    or if it's phased with the other parent
    four possible configurations of informative site and DNM that could
    provide phasing information...
    DNM---------------------INF
    alt---------------------alt >>> phased to parent with informative HET
    ref---------------------ref >>> phased to parent with informative HET
    alt---------------------ref >>> phased to other parent (spouse)
    ref---------------------alt >>> phased to other parent (spouse)
    """

    parent_evidence, spouse_evidence = 0, 0

    inf_ref_base, inf_alt_base = rd['inf_ref'], rd['inf_alt']
    dnm_ref_base, dnm_alt_base = rd['dnm_ref'], rd['dnm_alt']

    if inf_read_base == inf_alt_base and dnm_read_base == dnm_alt_base:
        return True
    elif inf_read_base == inf_ref_base and dnm_read_base == dnm_ref_base:
        return True
    elif inf_read_base == inf_alt_base and dnm_read_base == dnm_ref_base:
        return False
    elif inf_read_base == inf_ref_base and dnm_read_base == dnm_alt_base:
        return False

def phased(cram, chrom=None, inf_ref_pos=None, dnm_ref_pos=None, rd={}):
    """
    loop over every read aligned to the potential informative site,
    and for each read, catalog the base aligned to the informative
    and to the DNM site. if the position of the DNM lies outside of the
    read aligned to the informative site, check the mate of the read.
    """
    evidence = defaultdict(int)

    # loop over pileups in the BAM until we hit the position 
    # of the informative site to phase with
    read_cache = defaultdict()
    for pileupcolumn in cram.pileup(chrom, inf_ref_pos - 1, inf_ref_pos + 1):
        # do everything w/r/t the informative site
        if pileupcolumn.pos < inf_ref_pos: continue
        if pileupcolumn.pos > inf_ref_pos: break
        # for every read that spans the informative site, check that 
        # the read has support for both the informative HET and the known DNM.
        for pileupread in pileupcolumn.pileups:
            # ignore reads if there's an indel in the read, makes
            # it hard to adjust RBP. affects a small number of cases
            #if indel_in_read(pileupread.alignment) > 1: continue
            # ignore reads with indels w/r/t the reference
            if pileupread.is_del or pileupread.is_refskip: continue
            # skip reads that have MQ0
            inf_read_mq = pileupread.alignment.mapping_quality
            if inf_read_mq == 0: continue
            # get the aligned sequence
            inf_read_seq = pileupread.alignment.query_sequence
            # get the start and end of the aligned read w/r/t the reference
            # NOTE: these starts and ends do not include clipped bases
            inf_read_start = pileupread.alignment.reference_start
            inf_read_end = pileupread.alignment.reference_end
            # grab the position of the informative site in the query
            # sequence we're dealing with
            inf_read_pos = pileupread.query_position
            inf_read_base = inf_read_seq[inf_read_pos]
            # if the DNM isn't in the same read as the informative site,
            # we need to access the read's mate. however, accessing the mate
            # using pysam is inefficient, and moves us to a new position
            # in the CRAM file. so, instead, we'll cache the current read
            # (and the base aligned to the informative site) to use later
            if dnm_ref_pos not in range(inf_read_start, inf_read_end):
                # make sure the read has a mate that's aligned properly
                if not pileupread.alignment.is_proper_pair: continue
                if pileupread.alignment.query_name not in read_cache: 
                    read_cache[pileupread.alignment.query_name] = inf_read_base
                    continue
            else:
                # if the DNM *is* in the same read as the informative site,
                # get the position of the DNM in the read sequence
                dnm_read_base = get_dnm_base(cram, chrom, dnm_ref_pos, read=pileupread.alignment.query_name)
                if dnm_read_base is None: continue
                phased_to_inf_parent = phased_to_inf(rd, inf_read_base=inf_read_base, 
                                                         dnm_read_base=dnm_read_base)
                if phased_to_inf_parent:
                    evidence['inf'] += 1
                elif not phased_to_inf_parent:
                    evidence['other'] += 1

    if len(read_cache) == 0: return evidence

    # now, we loop back over every read aligned to the DNM position
    # and check if any of these reads are the mates of reads we cached
    # earlier. if so, we've cached the aligned base of that previous read.
    for pileupcolumn in cram.pileup(chrom, dnm_ref_pos - 1, dnm_ref_pos + 1):
        # only consider the DNM reference position
        if pileupcolumn.pos < dnm_ref_pos: continue
        if pileupcolumn.pos > dnm_ref_pos: break
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.query_name not in read_cache: continue
            if pileupread.is_del or pileupread.is_refskip: continue
            inf_read_base = read_cache[pileupread.alignment.query_name]
            dnm_read_mq = pileupread.alignment.mapping_quality
            if dnm_read_mq == 0: continue
            # NOTE: these starts and ends do not include clipped bases
            dnm_read_start, dnm_read_end = pileupread.alignment.reference_start, pileupread.alignment.reference_end
            # if the DNM isn't in the mate sequence, either, skip
            if dnm_ref_pos not in range(dnm_read_start, dnm_read_end): continue
            dnm_read_seq = pileupread.alignment.query_sequence
            # otherwise, grab the position of the DNM in the mate sequence
            dnm_read_pos = pileupread.query_position
            dnm_read_base = None 
            if pileupread.indel < 0:
                dnm_read_base = '-'
            elif pileupread.indel > 0:
                dnm_read_base = '+'
            else:
                try:
                    dnm_read_base = dnm_read_seq[dnm_read_pos]
                except IndexError: continue
            if dnm_read_base is None: continue
            phased_to_inf_parent = phased_to_inf(rd, inf_read_base=inf_read_base, 
                                          dnm_read_base=dnm_read_base)
            if phased_to_inf_parent:
                evidence['inf'] += 1
            elif not phased_to_inf_parent:
                evidence['other'] += 1

    return evidence

def get_evidence(dnms, samples, chrom=None,
                    prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):
    """
    loop over every one of the kid's variants in the BED file
    and get phasing evidence
    """
    import sys
    for i,d in enumerate(ts.reader(dnms, header="ordered")):
        kid = [s for s in samples if s.sample_id == d['sample_id']][0]
        # only loop over the specified range of DNMs
        if i % 100 == 0: print ("Done with {} variants".format(i), file=sys.stderr)
        if i == 0:
            keys = list(d.keys())
            keys.extend(['dad_evidence', 'mom_evidence'])
            print ('\t'.join(keys))
        if chrom and d['chrom'] != chrom: continue 
        # automatically phase male X DNMs to mom
        if d['chrom'] == 'X' and d['sample_sex'] == 'male':
            d['dad_evidence'] = '0'
            d['mom_evidence'] = '10'
            print ('\t'.join(d.values()))
            continue
        # get list of candidate parents and sites from the BED
        sites = d['candidate_sites'].split(',')
        if any([site == 'NA' for site in sites]): 
            d['dad_evidence'] = 'none'
            d['mom_evidence'] = 'none'
            print ('\t'.join(d.values()))
            continue
        combined_evidence = defaultdict(int)
        for site in sites:
            candidate_parent = site.split('|')[0]
            rest = site.split('|')[1]
            candidate_parent_sex = 'dad'
            if kid.mom.sample_id == candidate_parent: candidate_parent_sex = 'mom'
            inf_pos, dnm_pos = int(rest.split(':')[0]), int(d['start'])
            # get the reference and alternate base at the informative site
            inf_ref_base = rest.split(':')[-1].split('-')[0]
            inf_alt_base = rest.split(':')[-1].split('-')[1]
            # get the reference and alternate base at the DNM
            dnm_ref_base, dnm_alt_base = d['ref'], d['alt']
            # if we're dealing with a deletion, the 'pos' actually
            # needs to start at the end of the deletion
            del_len = None
            if len(dnm_ref_base) - len(dnm_alt_base) > 0:
                del_len = len(dnm_ref_base) - len(dnm_alt_base)
                dnm_alt_base = '-'
                dnm_ref_base = dnm_ref_base[0]
            elif len(dnm_alt_base) - len(dnm_ref_base) > 0:
                dnm_alt_base = '+'
            # create a dictionary of ref/alt bases for the informative site
            # and DNM (to pass into `phased`)
            rd = {'inf_ref':inf_ref_base, 'inf_alt':inf_alt_base,
                  'dnm_ref':dnm_ref_base, 'dnm_alt':dnm_alt_base, 
                  'del_len':del_len}
            # get read evidence for phasing from both parents
            bam = pysam.AlignmentFile(prefix + d['sample_id'] + '.bam', "rb")
            evidence = phased(bam, chrom=d['chrom'],
                                      inf_ref_pos=inf_pos,
                                      dnm_ref_pos=dnm_pos,
                                      rd=rd)
            if candidate_parent_sex == 'dad':
                combined_evidence['dad'] += evidence['inf']
                combined_evidence['mom'] += evidence['other']
            elif candidate_parent_sex == 'mom':
                combined_evidence['dad'] += evidence['other']
                combined_evidence['mom'] += evidence['inf']
        
        d['dad_evidence'] = str(combined_evidence['dad'])
        d['mom_evidence'] = str(combined_evidence['mom'])

        print ('\t'.join(d.values()))

def run(args):
    import os
    ped = Ped(args.ped)
    samples = [s for s in ped.samples()]
    print ("pysam version is {}".format(pysam.__version__), file=sys.stderr)
    dnms = gzip.open(args.dnms) if args.dnms.endswith('.gz') else args.dnms
    print ("starting to iterate through DNMs", file=sys.stderr)
    #fh = open(args.dnms + '.rbp', 'a')
    chrom = None
    if args.chrom:
        chrom = args.chrom
    get_evidence(dnms, samples, chrom=chrom)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    main(sys.argv[1:])

