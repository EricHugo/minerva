#! /usr/bin/env python

from __future__ import print_function
from collections import defaultdict, ChainMap
from Bio import SeqIO
from itertools import chain
from termcolor import cprint
import re

class findGeneNeighbourhood():
    def __init__(self, sequence, proteome, name, seq_length, matched_genes):
        self.name = name
        self.seq_length = seq_length
        self.sequence = sequence
        self.matched_genes = matched_genes
        with open(proteome) as protfile:
            self.pHeaders = set(header for header in protfile 
                    if re.search("^>", header))

    def evaluate_contig(self, match):
        found_n = False
        self.all_contig_loci = defaultdict(list)
        #print(match)
        for header in self.pHeaders:
            header = header.split(';')[0]
            #print(header)
            header_id = header.split('#')[-1].strip()
            #print(header_id)
            locus_tag = re.sub('>', '', header.split('#')[0].strip())
            contig_n = re.split("=|_", header_id)[1]
            #print(contig_n)
            self.all_contig_loci[contig_n].append(locus_tag)
            if match == locus_tag:
                found_n = contig_n
        if not found_n:
            cprint("something went horribly wrong", "red")
            cprint(self.name, "red")
            cprint(self.sequence, "red")
        #print(self.all_contig_loci[found_n])
        return self.all_contig_loci[found_n]
        
    def gather_locations(self, hmm):
        """Gathers all locations of all genes in a dict[gene] = [[START, STOP]]
        also populates a dict for the queried genes.
        """
        self.all_locs = defaultdict(list)
        self.hmm_locs = []
        all_contig_loci = self.evaluate_contig(hmm)
        for header in self.pHeaders:
            head = re.sub('>', '', header.split('#')[0].strip())
            #print(self.all_contig_loci)
            if head not in all_contig_loci:
                #print("continue")
                continue
            #print(head)
            self.all_locs[head].append(list(map(int, header.split('#')[1:3])))
            #print(hmm)
            if re.search(re.escape(hmm)+"\s", header):
                #cprint("match %s and %s" % (header, hmm), "red")
                self.hmm_locs.append(list(map(int, 
                    header.split('#')[1:3])))
        #print(self.all_locs)
        #print(self.hmm_locs)
        return self.hmm_locs, self.all_locs

    def check_overlap(self, flat_match, query_locs, reverse=False):
        """Tries to resolve an overlap in locations. Returns true if 
        starting query loc precedes the match start"""
        starts = flat_match[0] < query_locs[0]
        ends = flat_match[1] < query_locs[1]
        start_end = flat_match[0] > query_locs[1]
        end_start = flat_match[1] > query_locs[0]
        if reverse:
            return not bool(starts) and not bool(ends) and bool(end_start)
        return bool(starts) and bool(ends) and bool(end_start)

    def find_minimum_distance(self):
        """Compare the locations within hmm_locs to all_locs and find 
        gene with smallest positive and negative distance. Return 
        dict[hmm] = [minFloc, minRloc]
        """
        self.min_locs = defaultdict(list)
        for match, genes in self.matched_genes.items():
            match_locs, all_locs = self.gather_locations(match)
            #print(match_locs)
            #print(all_locs)
            flat_match = list(chain.from_iterable(match_locs))
            # try dict comp to include gene name with calculated distance
            # because else is mandatory in this syntax: None for values in
            # wrong direction, cannot simulate circularity since most genomes
            # are not closed
            # may implement method of determining whether to simulate circularity
            forward_locs = [{ head: int(loc[0] - flat_match[1]) 
                            if int(loc[0] - flat_match[1] > 0)
                            else 0 if self.check_overlap(flat_match, loc)
                            else None
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head == match ]
            # flatten list of dicts
            forward_dict = ChainMap(*forward_locs)
            reverse_locs = [{ head: int(flat_match[0] - loc[1]) 
                            if int(flat_match[0] - loc[1] > 0)
                            else 0 if self.check_overlap(flat_match, loc, True)
                            else None
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head == match ]
            reverse_dict = ChainMap(*reverse_locs)
            # min() is by far the slowest part of the code
            # inline lambda function actually faster than .get method
            #print(min(forward_dict, key=forward_dict.get))
            cprint(self.name, "yellow")
            try:
                # remove any None introduced neighbour detection
                forward_dict = { head: dist for head, dist in forward_dict.items() 
                                if type(dist) is int }
                forw_min = min(forward_dict.items(), key=lambda x: x[1])
                # some gbk have assembly gaps, that don't seperate into new 
                # contigs, therefore won't be caught as non-neighbour
                # therefore if gap is too long, assume its unreliable
                # and discard
                if forw_min[1] > 20000:
                    cprint(forward_locs, "red")
                    print(all_locs)
                    self.min_locs[match].append((None, None))
                else:
                    self.min_locs[match].append(forw_min)
            except ValueError:
                cprint(forward_dict, "cyan")
                self.min_locs[match].append((None, None))
            try: 
                reverse_dict = { head: dist for head, dist in reverse_dict.items() 
                                if type(dist) is int }
                rev_min = min(reverse_dict.items(), key=lambda x: x[1])
                if rev_min[1] > 20000:
                    cprint(reverse_locs, "red")
                    print(all_locs)
                    self.min_locs[match].append((None, None))
                else:
                    self.min_locs[match].append(rev_min)
            except ValueError:
                cprint(reverse_dict, "cyan")
                self.min_locs[match].append((None, None))
            cprint(self.min_locs, "green")
        return self.min_locs
