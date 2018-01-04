#! /usr/bin/env python

import pandas as pd
import multiprocessing as mp
import sys
import argparse
import numpy
import scipy
import csv

def _worker(task, DB, headers):
    # here we determine which function to call
    pass

def _listener(q):
    # we'll need a listener to safely output results from seperate
    # processes
    pass

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

    DB_all = pd.read_csv(args.minervaDB, sep='\t')

    # define tasks to be done, seperate each into multiprocessing
    queries = []
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads + 1)
    listener = pool.apply_async(_listener, (q,))

    jobs = []
    for i in queries:
        job = pool.apply_async(_worker, (i, DB_all, headers))
        jobs.append(job)

    for job in jobs:
        job.get()
    q.put("done")
    listener.get()
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
