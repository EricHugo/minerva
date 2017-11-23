#!/usr/bin/env python

from __future__ import print_function
from termcolor import cprint
import sys
import mmap
import re
import os


class parseMapFile():
    """The minerva database can be parsed using this module. Based on matching 
    regex patterns in specified column it will return the contents of the 
    requested column, by default it returns the path to the extracted gene 
    file.

    It is also capable of matching several queries in different columns. These 
    should be provided as lists with equal number of columns and queries.
    
    .. Example:  parseMapFile('minvera_db', ['match', 'CRISPR'], 
                              ['RAYT.*', 'True'], ['name'])

        This run would return the name of all genomes in minerva_db where any 
        RAYT is in the "match" column if the CRISPR column also contains 
        "True".
    """

    def __init__(self, map_file, query_column, queries, requests=["gene_path"], 
                 unique=False, unique_column=['name']):
        self.map_file = map_file
        self.query_column = query_column
        self.queries = queries
        self.unique = unique
        #print(self.queries)
        #print(self.query_column)
        if len(self.queries) is not len(self.query_column):
            raise AttributeError("Number of queries must be equal to number of \
                    columns")
        gene_files = open(self.map_file, 'r+b')
        self.gene_mem = mmap.mmap(gene_files.fileno(), 0, prot=mmap.PROT_READ)
        self.headers = str(self.gene_mem.readline(), "utf-8").split('\t')
        self.requests = requests
        self.col = [ col for col in self.find_column_num(self.query_column) ]
        self.req_col = [col for col in self.find_column_num(self.requests) ]
        if self.unique:
            self.uniq_cols = [ col for col in self.find_column_num(unique_column) ]
            self.unique_seen = set()
    
    def find_column_num(self, column_names):
        for name in column_names:
            i = 0
            for header in self.headers:
                if header.lower() == name.lower():
                    yield i
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
            if self.unique:
                uniq = self.check_uniq(mapping)
            else:
                uniq = True
            if not uniq:
                continue
            #print(unique_seen)
            pass_ = self.find_match(mapping)
            if pass_:
                request_list = [ mapping[r_col] for r_col in self.req_col ]
                return request_list
        self.gene_mem.seek(0)
        raise StopIteration

    def find_match(self, line):
        for sel, col in zip(self.queries, self.col):
            #print("match %s in %s" % (sel, col))
            if not re.fullmatch(sel, line[col], re.IGNORECASE):
                pass_ = False
                break
            #print("found %s in %s" % (sel, mapping[col]))
            pass_ = True
        return pass_

    def check_uniq(self, line):
        """Check that content of line in unique_columns have not previously
        been seen."""
        for uniq_col in self.uniq_cols:
            if line[uniq_col] in self.unique_seen:
                uniq = False
                break
            uniq = True
        # add to unique_seen only if found unique on all items (uniq=true)
        if uniq:
            for uniq_col in self.uniq_cols:
                self.unique_seen.add(line[uniq_col])
                #print(self.unique_seen)
        return uniq




