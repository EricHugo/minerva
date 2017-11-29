#! /usr/bin/env python

from __future__ import print_function
from collections import defaultdict, ChainMap
from itertools import chain
from termcolor import cprint
import re

class findGeneNeighbourhood():
    def __init__(self, proteome, name, seq_length, matched_genes):
        self.name = name
        self.seq_length = seq_length
        self.matched_genes = matched_genes
        with open(proteome) as protfile:
            self.pHeaders = set(header for header in protfile 
                    if re.search("^>", header))

    def gather_locations(self):
        """Gathers all locations of all genes in a dict[gene] = [[START, STOP]]
        also populates a dict for the queried genes.
        """
        self.all_locs = defaultdict(list)
        self.hmm_locs = defaultdict(list)
        for header in self.pHeaders:
            head = re.sub('>', '', header.split('#')[0].strip())
            #print(head)
            self.all_locs[head].append(list(map(int, header.split('#')[1:3])))
            for hmm, genes in self.matched_genes.items():
                #cprint(header, "green")
                #print(hmm)
                if re.search(re.escape(hmm)+"\s", header):
                    #cprint("match %s and %s" % (header, hmm), "red")
                    self.hmm_locs[hmm].append(list(map(int, 
                        header.split('#')[1:3])))
        #print(self.all_locs)
        #print(self.hmm_locs)
        return self.hmm_locs, self.all_locs

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
            #cprint(match_locs, "red")
            flat_match = list(chain.from_iterable(match_locs))
            # try dict comp to include gene name with calculated distance
            forward_locs = [{ head: int(loc[0] - flat_match[1]) 
                            if int(loc[0] - flat_match[1] > 0)
                            else int(loc[0] - flat_match[1] + self.seq_length) 
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head is match ]
            # flatten list of dicts
            forward_dict = ChainMap(*forward_locs)
            reverse_locs = [{ head: int(flat_match[0] - loc[1]) 
                            if int(flat_match[0] - loc[1] > 0)
                            else int(flat_match[0] - loc[1] + self.seq_length) 
                            for loc in locs } for head, locs in  
                            self.all_locs.items() if not head is match ]
            reverse_dict = ChainMap(*reverse_locs)
            # min() is by far the slowest part of the code
            # inline lambda function actually faster than .get method
            #print(min(forward_dict, key=forward_dict.get))
            forw_min = min(forward_dict.items(), key=lambda x: x[1])
            rev_min = min(reverse_dict.items(), key=lambda x: x[1])
            self.min_locs[match].append(forw_min)
            self.min_locs[match].append(rev_min)
            #print(self.min_locs)
        return self.min_locs
