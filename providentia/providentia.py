#! /usr/bin/env python

from contextlib import contextmanager
from termcolor import cprint
from scipy import stats
from collections import OrderedDict
import pandas as pd
import multiprocessing as mp
import multiprocessing.pool
import sys
import argparse
import numpy as np
import csv
import math
import time
import os
sys.path.insert(0, os.path.abspath(".."))

try:
    from recurse_database import recurse_database
except ImportError:
    from providentia.recurse_database import recurse_database


def _worker(df, init_header, headers, q, names_column, threads, depth):
    # here we determine how to recurse down the list
    # having received a init_header we start with significane testing it
    # then recurse down once level for each subgroup within init_header
    # -> testing -> recurse down -> testing down to len(header)
    # start looping back
    subgroups = df[init_header].unique()
    print(subgroups)
    print(threads)
    pool = mp.Pool(processes=threads)
    jobs = []
    for subgroup in subgroups:
        job = pool.apply_async(_init_run, (df, init_header, headers, subgroup,
                                           names_column, q))
        #job = pool.apply_async(__dummy, (subgroup, q))
        jobs.append(job)
        #print('started')
    for job in jobs:
        job.get()
    pool.close()
    pool.join()
    cprint("done?", "red")
    return 

def _init_run(df, init_header, headers, init_slice, names_column, q):
    """Initiates the recurse run, serves to allow multi-processing from 
    within the _worker child process."""
    sub_df = df[df[init_header].str.match(init_slice)]
    recurse_d = recurse_database(sub_df, init_header, headers, OrderedDict(
                                 {init_header: init_slice}), names_column, df,
                                 q=q)
    recurse_d.__next__()
    return

def __dummy(num, q):
    """Dummy process, repeatedly prints given variable and attempts to queue
    in given q"""
    while True:
        print(num)
        q.put(num)
        time.sleep(2)
        return

def _listener(q, outfile='-'):
    # we'll need a listener to safely output results from seperate
    # processes
    with open_stdout(outfile) as handle:
        while True:
            out_object = q.get()
            cprint(out_object, "green")
            if type(out_object) is tuple:
                stat = str(out_object[0])
                slices = str(out_object[1])
                handle.write(stat + ' ' + slices + '\n')
            try:
                if out_object == "done":
                    break
            except ValueError:
                # dataframes
                continue
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

def optimal_ints(total, div):
    """Produces an optimal list of integers given a total, a set size
    in descending order. If total < set size, populates remaining with
    1s."""
    opt_set = []
    while div: 
        quotient = total / div
        up_round = math.ceil(quotient)
        if up_round <= 0:
            # fill
            [ opt_set.append(1) for i in range(div) ]
            break
        div -= 1
        total -= up_round
        opt_set.append(up_round)
    return opt_set

class NoDaemonProcess(mp.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

class NoDaemonProcessPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess

def main():
    parser = argparse.ArgumentParser(description="""Attempts to find deviating 
            groups within a database produced based on all metrics present.
            Essentially performs an all v. all comparison.""")

    parser.add_argument("minervaDB", help="A database file produced by minerva")
    parser.add_argument("names_column", help="""The name of the column containing
                        unique identifiers per organism""")
    parser.add_argument("-c", "--columns", help="""Specify which columns will be
                        included in the analysis""")
    parser.add_argument("-e", "--entry_columns", help="""Specify which columns to
                        slice by initially.""")
    parser.add_argument("-d", "--depth", default=None, type=int,
                        help="""Specify maximum recursive depth that will be
                        descended""")
    parser.add_argument("-o", "--outfile", type=str, default='-', help="Specify \
                        outfile for results writing. By default writes to stdout.")
    parser.add_argument("-t", "--threads", default=1, type=int, help="""Define
                        number of threads to be used""")

    args = parser.parse_args()
    try:
        with open(args.columns) as cols:
            headers = cols.readline().strip().split('\t')
    except TypeError:
        with open(args.minervaDB) as DB:
            headers = DB.readline().strip().split('\t')
    print(headers)

    print(args.entry_columns)
    if args.entry_columns:
        with open(args.entry_columns) as entry:
            entry_set = entry.readline().strip().split('\t')
    df_minerva_all = pd.read_csv(args.minervaDB, sep='\t')
    names_column = args.names_column

    # define tasks to be done, seperate each into multiprocessing
    manager = mp.Manager()
    q = manager.Queue()
    # init pool as non-daemon to allow for spawning new child processes within
    pool = NoDaemonProcessPool(processes=args.threads + 1)
    listener = pool.apply_async(_listener, (q, args.outfile))

    jobs = []
    opt_threads = optimal_ints(args.threads, len(entry_set))
    for init_header, threads in zip(entry_set, opt_threads):
        job = pool.apply_async(_worker, (df_minerva_all, init_header, headers,
                                         q, names_column, threads, args.depth))
        jobs.append(job)
        #break

    for job in jobs:
        job.get()
    q.put("done")
    listener.get()
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
