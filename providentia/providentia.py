#! /usr/bin/env python

from contextlib import contextmanager
from scipy import stats
import pandas as pd
import multiprocessing as mp
import sys
import argparse
import numpy as np
import csv

def _worker(task, df_minerva, headers, selector, q):
    # here we determine which function to call
    # Outer_selection (e.g. RAYT1_pruned), inner_selection (e.g. Genus), target
    #  (e.g. Forward_distance)
    out = {}
    results = {}
    print(task)
    group, query = task
    print(group)
    # always remove any "-" entries as these represent no value in minerva
    print(len(df_minerva))
    df_valid = df_minerva[~df_minerva[query].str.match('-')].copy()
    print(len(df_valid))
    print(df_valid.dtypes)
    try:
        df_valid[query] = pd.to_numeric(df_valid[query], errors="raise")
    except ValueError:
        string_counts(df_valid, group, query)
        return
    df_valid.dropna(subset=[query], inplace=True)
    print(len(df_valid))
    print(df_valid.dtypes)
    all_values = df_valid[query].dropna().values
    #print(type(all_values))
    # here get all unique values in group column -> iterate over each
    # comparing with all_values in scipy ttest
    subgroups = df_valid[group].unique()
    for subgroup in subgroups:
        results[subgroup] = pandas_ttest(df_valid, all_values, group, query, 
                                        subgroup)
    print(results)
    out[selector] = results
    q.put(out)
    return

def string_counts(df, group, query_group):
    all_values = {}
    subgroups = df[group].unique()
    for subgroup in subgroups:
        sub_df = slice_dataframe(df, group, subgroup)
        print(subgroup)
        value = sub_df.groupby(query_group)[query_group].count()
        print(value.values)
        all_values[subgroup] = int(value.values)
    all_counts = all_values.values()
    print(all_counts)
    for subgroup, count in all_values.items():
        # can't simply find outlier counts, may have larger number of organisms
        # contributing to large counts
        # have to normalise
        pass
    return

def slice_dataframe(df, group, subgroup):
    sub_DB = df[df[group].str.match(subgroup)]
    return sub_DB

def transform_dataframe(df, query_group):
    # converts query column to numbers
    # either by counts or converting already numerical
    df_valid = df_minerva[~df_minerva[query_group].str.match('-')].copy()
    print(len(df_valid))
    print(df_valid.dtypes)
    try:
        df_valid[query_group] = pd.to_numeric(df_valid[query], errors="raise")
        df_valid.dropna(subset=[query_group], inplace=True)
        values = df_valid[query_group].dropna().values
    except ValueError:
        pass
    print(len(df_valid))
    print(df_valid.dtypes)
    return

def pandas_ttest(df, all_values, group, query, subgroup=None):
    if subgroup:
        sub_DB = df[df[group].str.match(subgroup)]
    else:
        sub_DB = df
    sub_values = sub_DB[query].values
    print(sub_values)
    ttest = stats.ttest_ind(all_values, sub_values, equal_var=False, 
                        nan_policy='raise')
    return ttest

def _listener(q, outfile='-'):
    # we'll need a listener to safely output results from seperate
    # processes
    with open_stdout(outfile) as handle:
        while True:
            out_object = q.get()
            if out_object == "done":
                break
            if out_object is type(dict):
                pass
    return

@contextmanager
def open_stdout(outfile):
    """Generator dynamically allows printing to file or stdout when used as open 
    with 'with'. No filename or '-' results in stdout, any other name will be
    used as filename to be written"""
    if outfile and outfile != '-':
        handle = open(outfile, 'w')
    else:
        handle = sys.stdout
    # close and flush open file if not stdout
    try:
        yield handle
    finally:
        if handle is not sys.stdout:
            handle.close()

def main():
    parser = argparse.ArgumentParser(description="""Attempts to find deviating 
            groups within a database produced based on all metrics present.
            Essentially performs an all v. all comparison.""")

    parser.add_argument("minervaDB", help="A database file produced by minerva")
    parser.add_argument("-t", "--threads", default=1, type=int, help="""Define 
                        number of threads to be used""")

    args = parser.parse_args()
    with open(args.minervaDB) as DB:
        headers = DB.readline().strip().split()
        print(headers)

    df_minerva_all = pd.read_csv(args.minervaDB, sep='\t')
    print(df_minerva_all.tail())
    outer_groups = ['RAYT1_pruned',]# 'RAYT2_pruned','RAYT3_pruned','RAYT4_pruned', 

    print(df_minerva_all.dtypes)

    # define tasks to be done, seperate each into multiprocessing
    tasks = [["Genus", "Match"]]
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads + 1)
    listener = pool.apply_async(_listener, (q,))

    jobs = []
    for outer_group in outer_groups:
        print(outer_group)
        if outer_group:
            outer_df_minerva = df_minerva_all.groupby('Match').get_group(outer_group)
        for task in tasks:
            job = pool.apply_async(_worker, (task, outer_df_minerva, headers, outer_group, q))
            jobs.append(job)

    for job in jobs:
        job.get()
    q.put("done")
    listener.get()
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
