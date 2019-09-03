#/usr/bin/env python

from scipy import stats
from collections import OrderedDict
from termcolor import cprint
from statsmodels.stats import proportion
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
            result = string_counts(parent_df, df, headers[-1], slices[headers[-1]])
            if result:
                result = (result, slices)
    return result

def num_rows_by_names(df, names, column, skip_nomatch=True):
    # TODO: Enable inclusion of no matches by skip_nomatch=False
    # TODO: Define "Match" column bar var
    ## here count "-" in column as 0
    ## but also, why?
    # try str.value_counts()
    copies = [len(df[df[column].str.match(name)]) if '-' not in
              df[df[column].str.match(name)]['Match'].tolist()
              else 0 for name in names]
    return copies

def subgroup(df, sub_df, column, query_group, names_column='Name'):
    all_groups = df[column].unique()
    if all_groups <= 1:
        return
    return

def string_counts(parent_df, df, group, query_group, names_column='Name',
                 kind=''):
    """Attempts to produce string counts per name in names_column
    compares the query_group to others in group to asses if it is
    outside 1 std dev."""
    print("string counts")
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
    analysis = copy_counts if kind == 'copy' else normalised_query_counts
    return analysis(df, subqueries, subgroups, names_column)

def copy_counts(df, queries, out_queries, names_column='Name'):
    """Function asks how many copies of Match by Name exists in query
    vs outside query. Significance returned if query median is +/-1 std
    from outside query median."""
    print("copy_counts")
    query_copies = num_rows_by_names(df, queries, names_column)
    group_copies = num_rows_by_names(df, out_queries, names_column)
    print(query_copies)
    print(group_copies)
    try:
        mean = np.mean(group_copies)
        std = np.std(group_copies)
        query_mean = np.mean(query_copies)
    except FloatingPointError:
        return
    print(mean, end=' ')
    print(std, end=' ')
    print(query_mean)
    if query_mean > mean + std or query_mean < mean - std:
        return (query_mean, mean, std)

def normalised_query_counts(df, query, out_queries, names_column='Name'):
    """Function asks how many matches are found under a query vs outside 
    queries normalised by total number of query and outside queries. Additional
    copies (matches under same names_column) are ignored."""
    # need to get all non-matches from query also
    print("norm")
    query_copies = [1 if not num == 0 else 0 for num in num_rows_by_names(df, query, names_column)]
    group_copies = [1 if not num == 0 else 0 for num in num_rows_by_names(df, out_queries, names_column)]
    print(query_copies)
    print(group_copies)
    try:
        q_prop = sum(query_copies) / len(query_copies)
        g_prop = sum(group_copies) / len(group_copies)
        stat, pval = proportion.proportions_ztest(q_prop, len(query_copies), g_prop)
    except (ZeroDivisionError, FloatingPointError) as E:
        cprint("norm error", "red")
        print(E)
        return
    print(stat, pval)
    if pval < 0.05:
        print(pval)
        return (q_prop, g_prop, len(query_copies))
    return

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
