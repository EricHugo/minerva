#/usr/bin/env python

from scipy import stats
from termcolor import cprint
from collections import OrderedDict
import pandas as pd
import numpy as np
import itertools
import sys
import random

try:
    from recurse_listener import recurse_listener
except ImportError:
    from .recurse_listener import recurse_listener

class recurse_database():
    def __init__(self, df, init_header, headers, init_slice, names_column, 
                 parent_df=None, max_depth=0, q=None):
        self.q = q
        try:
            if not parent_df:
                self.df = self.slice_df(df, init_header, init_slice)
                self.parent_df = df
        except ValueError:
            self.df = df
            self.parent_df = df
        self.init_header = init_header.split()
        self.headers = headers
        self.names = names_column
        self.listener = recurse_listener()
        self.slices = init_slice
        # initialise the listerner with next()
        next(self.listener)
        self.max_depth = max_depth
        self.current_depth = 0

    def __iter__(self):
        return self

    def __next__(self):
        return self.recurse(self.df, self.slices, self.init_header, self.parent_df)

    def __level__(self):
        return self.headerlist, self.slicelist

    def slice_df(self, df, column, query):
        sub_df = df[df[column].str.match(query)]
        return sub_df

    def get_headers(self, current_headers):
        """Using a list of headers, remove from total set and return a new 
        set with one more header"""
        num_headers = len(current_headers)
        #print(current_headers)
        headers = self.headers
        headers = [ head for head in headers if not head in current_headers ]
        #print("new")
        #print(headers)
        # currently uses permutations, which is a point of randomisation
        # avoid randomisation?
        header_permutations = [ header for header in itertools.permutations(headers, 1) ]
        #print(header_permutations)
        return header_permutations
        
    def recurse(self, df, slices, current_headers, parent_df=None):
        headers = self.get_headers(current_headers)
        #print("headers")
        #print(headers)
        ident = random.random()
        #print("new_")
        for header in headers:
            new_headers = current_headers + list(header)
            self.headerlist = new_headers
            subgroups = df[new_headers[-1]].unique()
            #print("headers")
            #print(new_headers)
            #print(ident)
            #print(slices)
            for subgroup in subgroups:
                slices[new_headers[-1]] = subgroup
                # here send to recurse_listener coroutine 
                # this will pause execution and handle results
                # resume upon completed yield loop
                try:
                    sub_df = self.slice_df(df, new_headers[-1], subgroup) 
                    df, df_type, result = self.call_listener(df, slices, new_headers,
                                                 parent_df)
                    if result:
                        self.q.put(result)
                except AttributeError:
                    # take a look at parent_df vs pre-slice df. should parent be
                    # the pre-sliced rather than from prev recurse
                    sub_df = df
                    df, df_type, result = self.call_listener(df, slices, new_headers,
                                                 parent_df)
                    # remove last slice. Necessary?
                    if result:
                            self.q.put(result)
                    del slices[new_headers[-1]]  
                    if df_type == "num":
                        break
                    continue
                # if there are not enough entries left in slice
                # or they are of numerical type
                # do not recurse
                if df_type == "num":
                    del slices[new_headers[-1]]
                    break
                    #continue
                #self.test.append(sub_df)
                self.recurse(df=sub_df, slices=slices, current_headers=new_headers,
                             parent_df=df)
                # clean last entry post-recurse
                del slices[new_headers[-1]]
        return 

    def call_listener(self,df, slices, new_headers, parent_df):
        """Calls the recurse_listener module"""
        valid_df, df_type = self.transform_dataframe(df, new_headers[-1])
        result = self.listener.send((valid_df, slices, new_headers, df_type, parent_df))
        return valid_df, df_type, result

    def transform_dataframe(self, df, query_group):
        """Removes empty values from query_group column ('-'), and attempts
        to convert the column to numerical. If not numerical returns type
        "str"."""
        try:
            # removes '-' entries, may not be desirable in all cases?
            ## add skip for when query is '-'
            # is necessary as long as we rely on on exception to determine
            # string / numeric
            try:
                df_valid = df[~df[query_group].str.match('-')].copy()
            except AttributeError:
                df_valid = df.copy()
            #print("received top")
            #print(query_group)
            #print(len(df_valid))
            #print(df_valid[query_group].dtypes)
            if df_valid[query_group].dtypes == "bool":
                #print("bool")
                raise ValueError
            # attempt to convert query column to numeric values
            df_valid[query_group] = pd.to_numeric(df_valid[query_group], 
                                                  errors="raise")
            # remove N/A from column
            df_valid.dropna(subset=[query_group], inplace=True)
            # check if actual values remain in the query group
            values = df_valid[query_group].dropna().values
            # if there are no values, implies dataframe for the group is only '-'
            if not values.any():
                #print("no values")
                #print(len(df))
                #print(df[query_group])
                # if we raise ValueError here -> will be treated as strings
                # therefore num of copies
                raise ValueError
            #print(values)
            df_type = "num"
        except ValueError:
            df_valid = df.copy()
            df_valid[query_group] = df[query_group].astype(str)
            df_type = "str"
            #print("received exception")
        #print(len(df_valid))
        #print(df_valid.dtypes)
        #print(df_valid)
        return df_valid, df_type


