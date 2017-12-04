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
        """Checks which contig any given match is on, and extracts only 
        other CDS on the same contig. This prevents erroneous distance
        measures"""
        found = False
        #print(match)
        with open(self.sequence) as handle:
            for record in SeqIO.parse(handle, 'genbank'):
                for feature in record.features:
                    if feature.type == 'source':
                        if found:
                            #print(self.all_contig_loci)
                            return self.all_contig_loci
                        self.all_contig_loci = []
                    if feature.type == 'CDS':
                        self.all_contig_loci.append(''.join(
                                feature.qualifiers['locus_tag']))
                        if match == self.all_contig_loci[-1]:
                            found = True
        if not found:
            cprint("something went horribly wrong", "red")
            cprint(self.name, "red")
            cprint(self.sequence, "red")
        return None

    def gather_locations(self):
        """Gathers all locations of all genes in a dict[gene] = [[START, STOP]]
        also populates a dict for the queried genes.
        """
        self.all_locs = defaultdict(list)
        self.hmm_locs = defaultdict(list)
        for hmm, genes in self.matched_genes.items():
            self.evaluate_contig(hmm)
            for header in self.pHeaders:
                head = re.sub('>', '', header.split('#')[0].strip())
                #print(self.all_contig_loci)
                if head not in self.all_contig_loci:
                    #print("continue")
                    continue
                #print(head)
                self.all_locs[head].append(list(map(int, header.split('#')[1:3])))
                #print(hmm)
                if re.search(re.escape(hmm)+"\s", header):
                    #cprint("match %s and %s" % (header, hmm), "red")
                    self.hmm_locs[hmm].append(list(map(int, 
                        header.split('#')[1:3])))
        ##print(self.all_locs)
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
        try:
            self.hmm_locs
        except AttributeError:
            self.gather_locations()
        self.min_locs = defaultdict(list)
        for match, match_locs in self.hmm_locs.items():
            flat_match = list(chain.from_iterable(match_locs))
            # try dict comp to include gene name with calculated distance
            forward_locs = [{ head: int(loc[0] - flat_match[1]) 
                            if int(loc[0] - flat_match[1] > 0)
                            else 0 if self.check_overlap(flat_match, loc)
                            else int(loc[0] - flat_match[1] + self.seq_length) 
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head == match ]
            # flatten list of dicts
            forward_dict = ChainMap(*forward_locs)
            reverse_locs = [{ head: int(flat_match[0] - loc[1]) 
                            if int(flat_match[0] - loc[1] > 0)
                            else 0 if self.check_overlap(flat_match, loc, True)
                            else int(flat_match[0] - loc[1] + self.seq_length) 
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head == match ]
            reverse_dict = ChainMap(*reverse_locs)
            # min() is by far the slowest part of the code
            # inline lambda function actually faster than .get method
            #print(min(forward_dict, key=forward_dict.get))
            forw_min = min(forward_dict.items(), key=lambda x: x[1])
            rev_min = min(reverse_dict.items(), key=lambda x: x[1])
            self.min_locs[match].append(forw_min)
            self.min_locs[match].append(rev_min)
            cprint(self.min_locs, "green")
        return self.min_locs
