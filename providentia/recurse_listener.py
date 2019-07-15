#/usr/bin/env python

from scipy import stats
from collections import OrderedDict
from termcolor import cprint
import sys
import pandas as pd
import numpy as np
import time
import warnings

np.seterr(all='raise')

def recurse_listener():
    """Coroutine to recurse_database.py, listens to sent objects containing
    a DataFrame and slice information."""
    result = None
    while True:
        df, slices, headers, df_type, parent_df = yield result
        #print("df_type: " + df_type)
        # names column needs to be sent also, determine number of copies
        names = None
        # here call appropriate stats test
        # send results to queue(?) => listener
        #print(slices)
        #print(headers)
        if df_type == "num":
            try:
                result = pandas_ttest(parent_df, headers[-1], slices[headers[-2]],
                                      headers[-2])
                if result:
                    result = (result, slices)
            except KeyError:
                print("KeyError")
                print(slices)
                print(headers)
                continue
        elif df_type == "str":
            #print("str")
            print(headers)
            print(slices)
            result = copy_numbers(parent_df, df, headers[-1], slices[headers[-1]])
            if result:
                result = (result, slices)
    return result

def copy_numbers(parent_df, df, group, query_group, names_column='Name'):
    """Attempts to produce string counts per name in names_column
    compares the query_group to others in group to asses if it is
    outside 1 std dev."""
    print("copy counts")
    print(df)
    print(group)
    print(query_group)
    try:
        subgroups = df[names_column].unique()
    except TypeError:
        # not sure why typeerror. Empty df?
        cprint("copy_numbers() TypeError", "red")
        print(df)
        return
    # check if there are enough values to do any population comps
    # else skip
    if len(subgroups) <= 1:
        return
    sub_df = slice_dataframe(df, group, query_group)
    #print("sub")
    print(sub_df)
    subqueries = sub_df[names_column].unique()
    subgroups = [subgroup for subgroup in subgroups
                 if subgroup not in subqueries]
    print("subgroups: ")
    print(subgroups)
    ## here count "-" in column as 0
    ## but also, why?
    query_copies = [len(df[df[names_column].str.match(name)]) if '-' not in
                    df[df[names_column].str.match(name)]['Match'].tolist()
                    else 0 for name in subqueries]
    group_copies = [len(df[df[names_column].str.match(name)]) if '-' not in
                    df[df[names_column].str.match(name)]['Match'].tolist()
                    else 0 for name in subgroups]
    print(query_copies)
    print(group_copies)
    try:
        mean = np.mean(group_copies + query_copies)
        std = np.std(group_copies + query_copies)
        query_mean = np.mean(query_copies)
    except FloatingPointError:
        return
    print(mean, end=' ')
    print(std, end=' ')
    print(query_mean)
    if query_mean > mean + std or query_mean < mean - std:
        return (query_mean, mean, std)

def slice_dataframe(df, group, subgroup):
    """Takes a dataframe, group and subgroup of group.
    Returns dataframe with only entries matching subgroup"""
    sub_DB = df[df[group].str.match(subgroup)]
    #print(sub_DB)
    return sub_DB

def df_to_num(df, column):
    try:
        df_valid = df[~df[column].str.match('-')].copy()
    except AttributeError:
        df_valid = df.copy()
    #print("received top")
    #print(column)
    #print(len(df_valid))
    #print(df_valid[column].dtypes)
    if df_valid[column].dtypes == "bool":
       #print("bool")
        raise ValueError
    # attempt to convert query column to numeric values
    df_valid[column] = pd.to_numeric(df_valid[column], errors="raise")
    # remove N/A from column
    df_valid.dropna(subset=[column], inplace=True)
    #print(df_valid)
    return df_valid

def pandas_ttest(parent_df, query, target, target_group, names_column="Name",
                 cutoff=0.1):
    """Execute a simple ttest for values in query column as compared to other
    valeus in column"""
    # correct for duplicate Match entries here?
    # otherwise shifts population due to num matches, not genomes
    df = df_to_num(parent_df, query)
    sub_df = slice_dataframe(df, target_group, target)
    sub_values = sub_df[query].values
    # anti-match the sub_df from wider df to seperate all from sub
    merged_df = pd.merge(df, sub_df, indicator=True, how='outer')
    outer_df = merged_df[merged_df['_merge'].str.match('left_only')]
    all_values = outer_df[query].values
    #print("ttest here")
    #print("query: " + query)
    #print(all_values)
    #print(sub_values)
    #print(target)
    #print(target_group)
    #print(df)
    # preform Welch's t-test
    # potential problem with 0-distance values
    if not sub_values and all_values:
        return
    try:
        ttest = stats.ttest_ind(all_values, sub_values, equal_var=False, 
                                nan_policy='raise')
    except FloatingPointError:
        cprint("warning for:", "magenta")
        print(all_values)
        print(sub_values)
    #print(ttest)
    if ttest.pvalue <= cutoff:
        return ttest
    #print("query: " + query, end=' ')
    #print("target: " + target)
    #print(ttest)
    #sys.exit()
