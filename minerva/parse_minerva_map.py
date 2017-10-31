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
        gene_files = open(self.map_file, 'r+b')
        self.gene_mem = mmap.mmap(gene_files.fileno(), 0, prot=mmap.PROT_READ)
        headers = str(self.gene_mem.readline(), "utf-8").split('\t')
        self.request = request
        self.col = 0
        self.req_col = 0
        for header in headers:
            if header.lower() == self.query_column.lower():
                break
            self.col += 1
        for header in headers:
            if header.lower() == self.request.lower():
                break
            self.req_col += 1

    def __iter__(self):
        return self

    def __next__(self):
        return self.gene_paths()

    def gene_paths(self):
        """Extracts and outputs iterable of all gene pathes matching query in 
        column"""
        for mapping in iter(self.gene_mem.readline, bytes()):
            mapping = str(mapping, "utf-8").split('\t')
            if self.query_selection in mapping[self.col]:
                return mapping[self.req_col]
        self.gene_mem.seek(0)
        raise StopIteration


