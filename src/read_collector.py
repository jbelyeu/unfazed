#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import pysam
import numpy as np
MILLION=1000000
MIN_MAPQ=1
STDEV_COUNT=3

def estimate_discordant_insert_len(bamfile):
    insert_sizes = []
    for i,read in enumerate(bamfile):
        insert_sizes.append(abs(read.tlen))
        if i >= MILLION:
            break
    insert_sizes = np.array(insert_sizes)
    #filter out the top .1% as weirdies
    insert_sizes = np.percentile(insert_sizes, 99.5)
    
    frag_len = int(np.mean(insert_sizes))
    stdev = np.std(insert_sizes)
    return frag_len+(stdev*STDEV_COUNT)

def goodread(read):
    if (read.is_qcfail
            or read.is_unmapped
            or read.is_duplicate
            or int(read.mapping_quality) < MIN_MAPQ
            or read.is_secondary
            or read.is_supplementary
            or read.mate_is_unmapped
            or (read.next_reference_id != read.reference_id)
    ):
        return False
    return True


def collect_reads_sv(bam_name, region, discordant_len=None):
    """
    given an alignment file name, and a de novo SV region,
    return the reads that support the variant as a dictionary with two lists,
    containing reads that support and reads that don't (latter will be none in an SV)
    """
    bamfile = pysam.AlignmentFile(bam_name, 'rb')

    if not discordant_len:
        discordant_len = estimate_discordant_insert_len(bamfile)

    supporting_reads = []
    for position in region['start'],region['end']:
        position = int(position)
        bam_iter = bamfile.fetch(
            region['chrom'], 
            position-discordant_len, 
            position+discordant_len
        )
        readcount = 0
        for read in bam_iter:
            if not goodread(read):
                continue
            if (read.has_tag("SA")):
                #if it's a splitter and the softclipping begins within 5 bases of the breaks, keep it
                if (((position == int(region['start'])) and ((position-5) <= read.get_reference_positions()[-1] <= (position+5))) 
                      or ((position == int(region['end'])) and ((position-5) <= read.get_reference_positions()[0] <= (position+5)))):
                    supporting_reads.append(read)

                    #find mate for informative site check
                    mate = bamfile.mate(read)

                    if goodread(mate):
                        supporting_reads.append(mate)
    
            elif read.tlen > discordant_len:
                #find mate for informative site check
                mate = bamfile.mate(read)
                read_positions = [
                    mate.reference_start,
                    read.reference_start,
                    mate.reference_end,
                    read.reference_end
                ]
                if (goodread(mate) and min(read_positions) <= int(region['start']) <= max(read_positions)):
                    supporting_reads.append(mate)
                    supporting_reads.append(read)
    return {
        "alt" : supporting_reads,
        "ref" : []
    }

if __name__ == "__main__":
    sys.exit("Import this as a module")