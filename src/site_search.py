#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys

def binary_search(start, end, informative_sites_list):
    match = None
    query_start = 0
    query_end = len(informative_sites_list)

    while (not match and query_end > 0):
        #check the middle
        query_pos = int((query_end-query_start)/2)
        if (query_start - query_end) == 0:
            break

        #if the query position is between the start and the end of the read, we've arrived
        if start <= informative_sites_list[query_pos]['pos'] <= end:
            match = informative_sites_list[query_pos]
            break
        elif informative_sites_list[query_pos]['pos'] > start: 
            #move left, position is too high
            query_end = query_pos-1
        elif informative_sites_list[query_pos]['pos'] < start: 
            #move right, position is too low
            query_start = query_pos+1
        else:
            print("this should be impossible")
    return match


def match_informative_sites(reads, informative_sites, chrom):
    """
    Given a list of pysam reads, 
    and a dictionary of lists of informative sites, 
    return the informative sites that are in the reads.
    Pasing in chrom because the pysam object keeps getting it wrong
    """
    matches = {}

    #reads are stored by whether they support ref or alt alleles
    for ref_alt in reads:
        matches[ref_alt] = []

        for read in reads[ref_alt]:
            site_match = binary_search(
                read.reference_start, 
                read.reference_end, 
                informative_sites[chrom]
            )
            if site_match:
                site_match['read'] = read
                matches[ref_alt].append(site_match)
    return matches

if __name__ == "__main__":
    sys.exit("Import this as a module")
