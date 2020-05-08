#! /usr/bin/env python
from __future__ import print_function

# Python 2/3 compatibility
import sys


def binary_search(start, end, informative_sites):
    matches = []
    query_start = 0
    query_end = len(informative_sites) - 1
    query_start_prev = -1
    query_end_prev = -1
    while len(matches) <= 0 and query_end > -1:
        # if the query region converges, no sites in read
        if query_start > query_end:
            break
        # if the query region isn't changing, no sites in read
        if (query_start == query_start_prev) and (query_end == query_end_prev):
            break

        query_start_prev = query_start
        query_end_prev = query_end
        query_pos = int((query_end + query_start) / 2)
        # if the query position is between the start and the end of the read, we've arrived
        # print(query_start, query_pos, query_end)
        if start <= informative_sites[query_pos]["pos"] < end:
            # keep the one that we found
            matches.append(informative_sites[query_pos])

            # check the next ones to see if they also fit the bill
            for isite in informative_sites[query_pos + 1:]:
                if start <= isite["pos"] <= end:
                    matches.append(isite)
                else:
                    break
            # and the previous ones
            for isite in informative_sites[:query_pos][::-1]:
                if start <= isite["pos"] <= end:
                    matches.append(isite)
                else:
                    break
            break
        elif informative_sites[query_pos]["pos"] > start:
            # move left, position is too high
            query_end = query_pos - 1
        elif informative_sites[query_pos]["pos"] < start:
            # move right, position is too low
            query_start = query_pos + 1
    return matches


def match_informative_sites(reads, informative_sites):
    """
    Given a list of pysam reads,
    and a dictionary of lists of informative sites,
    return the informative sites that are in the reads.
    Pasing in chrom because the pysam object keeps getting it wrong
    """
    matches = {}

    # reads are stored by whether they support ref or alt alleles
    for ref_alt in reads:
        matches[ref_alt] = []

        for read in reads[ref_alt]:
            site_matches = binary_search(
                read.reference_start, read.reference_end, informative_sites
            )
            if len(site_matches) > 0:
                match_info = {"matches": site_matches, "read": read}
                matches[ref_alt].append(match_info)
    return matches


if __name__ == "__main__":
    sys.exit("Import this as a module")
