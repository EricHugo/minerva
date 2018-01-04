#! /usr/bin/env python

import pandas as pd
import sys
import argparse
import numpy
import scipy
import csv

def _worker(DB, headers):
    # here we define which 
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
        headers = DB.readline()
        print(headers)

    DB_all = pd.read_csv(args.minervaDB, sep='\t')

if __name__ == "__main__":
    main()
