#/usr/bin/env python

from scipy import stats
import pandas as pd
import numpy as np
import itertools
import sys

class recurse_database():
    def __init__(self, df, init_header, headers, names_column):
        self.df = df
        self.init_header = init_header.split()
        self.headers = headers
        self.names = names_column

    def __iter__(self):
        return self

    def __next__(self):
        return self.recurse(self.df, self.init_header)

    def __level__(self):
        return self.headerlist, self.slicelist

    def slice_df(self, df, column, query):
        sub_df = df[df[column].str.match(query)]
        return sub_df

    def get_headers(self, current_headers):
        num_headers = len(current_headers)
        print(current_headers)
        headers = self.headers
        headers = [ head for head in headers if not head in current_headers ]
        for header_permutation in itertools.permutations(headers, 1):
            yield header_permutation 
        return

    def recurse(self, df, current_headers):
        for header in self.get_headers(current_headers):
            new_headers = current_headers + list(header)
            print(new_headers)
            print(header)
            self.headerlist = new_headers
            subgroups = df[new_headers[-1]].unique()
            print(subgroups)
            for subgroup in subgroups:
                self.slicelist = None
                sub_df = self.slice_df(df, new_headers[-1], subgroup) 
                self.recurse(df=sub_df, current_headers=new_headers)
                return sub_df
