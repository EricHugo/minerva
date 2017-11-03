#!/usr/bin/env python

from __future__ import print_function
import sys
import mmap
import re
import os


class parseMapFile():
    def __init__(self, map_file, query_column, query_selection, request="gene_path"):
        self.map_file = map_file
        self.query_column = query_column
        self.query_selection = query_selection
        #print(self.query_selection)
        #print(self.query_column)
        if len(self.query_selection) is not len(self.query_column):
            raise AttributeError("Number of queries must be equal to number of \
                    columns")
        gene_files = open(self.map_file, 'r+b')
        self.gene_mem = mmap.mmap(gene_files.fileno(), 0, prot=mmap.PROT_READ)
        headers = str(self.gene_mem.readline(), "utf-8").split('\t')
        self.request = request
        self.col = []
        self.req_col = []
        for q in self.query_column:
            i = 0
            for header in headers:
                if header.lower() == q.lower():
                    self.col.append(i)
                    break
                i += 1
        for r in self.request:
            i = 0
            for header in headers:
                if header.lower() == r.lower():
                    self.req_col.append(i)
                    break
                i += 1

    def __iter__(self):
        return self

    def __next__(self):
        return self.gene_paths()

    def gene_paths(self):
        """Extracts and outputs iterable of all gene pathes matching query in 
        column"""
        for mapping in iter(self.gene_mem.readline, bytes()):
            mapping = str(mapping, "utf-8").split('\t')
            for sel, col in zip(self.query_selection, self.col):
                #print("match %s in %s" % (sel, col))
                if not re.fullmatch(sel, mapping[col], re.IGNORECASE):
                    pass_ = False
                    break
                #print("found %s in %s" % (sel, mapping[col]))
                pass_ = True
            if pass_:
                request_list = [ mapping[r_col] for r_col in self.req_col ]
                return request_list
        self.gene_mem.seek(0)
        raise StopIteration


