#!/usr/bin/env python

from __future__ import print_function
from matplotlib import pyplot as plt
from collections import defaultdict
import matplotlib
import numpy as np
import argparse
import sys
import os
import copy

sys.path.append('/home/hugoson/git/minerva')
try:
    from minerva import parseTaxonomy, parseMapFile
except ImportError:
    from parse_taxonomy import parseTaxonomy
    from parse_minerva_map import parseMapFile

matplotlib.style.use('ggplot')

def _worker(minerva_table, search_set):
    found_dict = defaultdict(list)
    search_lists = [ search.split(' ') for search in search_set ]
    print(search_lists)
    found_list = [ ''.join(found) for found in parseMapFile(minerva_table, 
        search_lists[0], search_lists[1], search_lists[2]) ]
    #print(found_list)
    all_search_lists = copy.deepcopy(search_lists)
    all_search_lists[0][0] = "name"
    all_search_lists[1][0] = ".*"
    all_found_list = [ ''.join(all_found) for all_found in parseMapFile(minerva_table, 
        all_search_lists[0], all_search_lists[1], all_search_lists[2]) ]
    #print(all_found_list)
    found_dict[search_lists[1][1]].append(found_list)
    found_dict[search_lists[1][1]].append(all_found_list)
    try:
        secondary_selection = search_lists[1][0]
    except IndexError:
        secondary_selection = None
    #print(found_dict)
    return found_dict, secondary_selection

def create_bar_summary(found_v_total_set, secondary_title):
    fig, ax = plt.subplots()
    sep = 0
    for found_v_total in found_v_total_set:
        #print(found_v_total)
        found_no = [ len(set(found[0])) for name, found in found_v_total.items() ]
        remaining_no = [ len(set(total[1])) - len(set(total[0])) for name, total in 
                found_v_total.items() ]
        labels = [ name for name, found in found_v_total.items() ]
        ind = np.arange(len(found_v_total))
        p1 = ax.bar(ind + sep, remaining_no, color='#348abd', width=0.15, )
        p2 = ax.bar(ind + sep, found_no, width=0.15, color='#e24a33', 
                bottom=remaining_no)
        sep += 0.19
    ax.set_xticks(ind + sep / 2)
    ax.set_xticklabels(labels, rotation=45, fontsize=9, ha='right')
    plt.ylabel('Taxa')
    #
    plt.legend((p1[0], p2[0]), ('Total', 'Found'))
    plt.show()
    return

def main():
    parser = argparse.ArgumentParser(description="""Takes result table from
    minerva and summarises results by desired slice.""")

    parser.add_argument("minervaTable", help="""Results table produced by
    minerva""")
    parser.add_argument("-q", "--query", required=False, nargs='*', help="""Query 
    or queries (space seperated) to be searched for in columns specified""")
    parser.add_argument("-c", "--column", required=False, nargs='+', 
    help="""Columns to be searched through for specified queries. Order sensitive""")
    parser.add_argument("-r", "--request", required=False, nargs='*', help="""
    Column of matched to output""")
    parser.add_argument("-f", "--file", required=False, help="""File in which 
    mutiple seperate queries and columns summaries are requested. Structured with 
    seperate queries with space, queries and columns seperated by tab, and 
    seperated searches by newline.""")

    args = parser.parse_args()
    
    if args.file:
        with open(args.file) as queries:
            found_series = defaultdict(list)
            found_series_set = []
            old_ss = None
            for each in queries:
                found_set, secondary_selection = _worker(args.minervaTable, 
                        each.strip().split('\t'))
                print(secondary_selection)
                if secondary_selection != old_ss and old_ss:
                    found_series_set.append(found_series)
                    found_series = defaultdict(list)
                found_series = {**found_series, **found_set}
                old_ss = secondary_selection
            found_series_set.append(found_series)
            print(len(found_series_set))
    create_bar_summary(found_series_set, secondary_selection)

if __name__ == "__main__":
    main()
