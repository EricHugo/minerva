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
    group, query = task
    print(group)
    df_minerva[query] = pd.to_numeric(df_minerva[query], errors="coerce")
    df_minerva.dropna(subset=[query], inplace=True)
    #print(df_minerva.dtypes)
    all_values = df_minerva[query].dropna().values
    #print(type(all_values))
    # here get all unique values in group column -> iterate over each
    # comparing with all_values in scipy ttest
    subgroups = df_minerva[group].unique()
    for subgroup in subgroups:
        results[subgroup] = pandas_ttest(df_minerva, group, subgroup)
    print(results)
    out[selector] = results
    q.put(out)
    return

def pandas_ttest(df, group, subgroup=None):
    if subgroup:
        sub_DB = df[df[group].str.match(subgroup)]
    else:
        sub_DB = df
    sub_values = sub_DB[query].values
    print(sub_values)
    t = stats.ttest_ind(all_values, sub_values, equal_var=False, 
                        nan_policy='raise')
    print(t)
    return t


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
    outer_groups = ['RAYT1_pruned']#, 'RAYT2_pruned','RAYT3_pruned','RAYT4_pruned']

    # define tasks to be done, seperate each into multiprocessing
    tasks = [["Genus", "Forward_distance"]]
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads + 1)
    listener = pool.apply_async(_listener, (q,))

    jobs = []
    for outer_group in outer_groups:
        print(outer_group)
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
