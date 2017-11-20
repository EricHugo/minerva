#!/usr/bin/env python

from __future__ import print_function
from matplotlib import pyplot as plt
from collections import defaultdict, OrderedDict
from itertools import chain
from termcolor import cprint
import pandas as pd
import seaborn as sns
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

def gather_minerva_data(minerva_table, search_lists, gather_negatives=False):
    all_found_list = []
    selectors = []
    print(search_lists)
    # extract list of matches based on all given columns and queries per given 
    # column
    for search_list in search_lists:
        found_dict = defaultdict(list)
        cprint(search_list, "red")
        query_list = [ query.split(' ') for query in search_list ]
        print(query_list[2])
        found_list = [ found for found in parseMapFile(minerva_table, 
            query_list[0], query_list[1], query_list[2], unique=True, 
            ) ]
        print(found_list)
        found_dict[query_list[1][-1]].append(found_list)
        all_found_list.append(found_dict)
        selectors.append(query_list[1][0])
        print(selectors)
    # extract all results, assuming it matches index 1+, to get a total possible
    #print(all_found_list)
    # specify last index of index 1 as the name for the particular list
    # find the primary selector of the query
    # if there are more selectors, grab the last selector and return
    if len(query_list[0]) > 1:
        secondary_selector = query_list[0][-1]
    else:
        secondary_selector = None
    #print(found_dict)
    print(secondary_selector)
    return all_found_list, selectors, secondary_selector

def _worker(minerva_table, search_set):
    found_dict = defaultdict(list)
    # split each section by space, treat as queries
    search_lists = [ search.split(' ') for search in search_set ]
    #print(search_lists)
    # extract list of matches based on all given columns and queries per given 
    # column
    found_list = [ ''.join(found) for found in parseMapFile(minerva_table, 
        search_lists[0], search_lists[1], search_lists[2]) ]
    #print(found_list)
    # deepcopy the query list so it may be altered
    all_search_lists = copy.deepcopy(search_lists)
    # replace first column and query with matching every single name
    all_search_lists[0][0] = "name"
    all_search_lists[1][0] = ".*"
    # extract all results, assuming it matches index 1+, to get a total possible
    all_found_list = [ ''.join(all_found) for all_found in parseMapFile(minerva_table, 
        all_search_lists[0], all_search_lists[1], all_search_lists[2]) ]
    #print(all_found_list)
    # specify last index of index 1 as the name for the particular list
    found_dict[search_lists[1][-1]].append(found_list)
    found_dict[search_lists[1][-1]].append(all_found_list)
    # find the primary selector of the query
    selector = search_lists[1][0]
    # if there are more selectors, grab the last selector and return
    try:
        secondary_selector = search_lists[0][-1]
    except IndexError:
        secondary_selector = None
    #print(found_dict)
    return found_dict, selector, secondary_selector

def create_stacked_summary(found_v_total_set, secondary_title):
    fig, ax = plt.subplots()
    sep = 0
    width = 0.7 / len(found_v_total_set)
    print(type(found_v_total_set))
    for each in found_v_total_set:
        print([ ind for ind, val in each.items() ])
    for found_v_total in found_v_total_set:
        #print(found_v_total)
        found_no = [ len(set(found[0])) / len(set(found[1])) for name, found in found_v_total.items() ]
        remaining_no = [ (len(set(total[1])) - len(set(total[0]))) / len(set(total[1])) for name, total in 
                found_v_total.items() ]
        labels = [ name for name, found in found_v_total.items() ]
        ind = np.arange(len(found_v_total))
        print(labels)
        print(ind)
        p1 = ax.bar(ind + sep, remaining_no, color='#348abd', width=width)
        p2 = ax.bar(ind + sep, found_no, width=width, color='#e24a33', 
                bottom=remaining_no)
        sep += width + 0.05
    # set xticks to be middle for each set, subtracts width and 0.05 to
    # compensate last increment
    ax.set_xticks(ind + (sep - width - 0.05) / 2)
    ax.set_xticklabels(labels, rotation=45, fontsize=9, ha='right')
    plt.ylabel('Taxa')
    #
    plt.legend((p1[0], p2[0]), ('Negative', 'Postive'))
    plt.show()
    return

def create_frequency_summary(minerva_table, search_sets):
    # try pandas and seaborn for grouping
    df = pd.DataFrame()
    all_found_matches, selectors, secondary_selector = gather_minerva_data(
            minerva_table, search_sets)
    for selector, found_v_total in zip(selectors, all_found_matches):
        print(selector)
        data = {}
        # questionable whether the if clause should be included, possibly better 
        # to handle empty list counts
        found_all = [ found[0] for name, found in found_v_total.items() if found[0] ]
        labels =  [ name for name, found in found_v_total.items() if found[0]] 
        #print(found_all)
        # necessary because nested list
        for found, label in zip(found_all, labels):
            print(label)
            print(found)
            found_flat = list(chain(*found))
            try:
                found_numbers = [ float(num) for num in found_flat ]
            except ValueError:
                found_numbers = [ found_flat.count(num) for num in set(found_flat) ] 
            if secondary_selector:
                second_sel = [ label for i in range(len(found_numbers)) ]
                data[secondary_selector] = second_sel
            prim_sel = [ selector for i in range(len(found_numbers)) ]
            data["values"] = found_numbers
            data["Selector"] = prim_sel
            #print(found_numbers) 
            #print(label)
            df_new = pd.DataFrame(data)
            df = df.append(df_new, ignore_index=True)
    print(df)
    print(secondary_selector)
    if secondary_selector:
        bp = sns.boxplot(x='Selector',y='values',data=df,hue=secondary_selector, width=0.7)
    else:
        bp = sns.boxplot(x='Selector',y='values',data=df, width=0.7)
    plt.show()
    return

def create_residual_summary(minerva_table, search_sets):
    df = pd.DataFrame()
    search_lists = []
    print(search_sets)
    # 
    all_found_matches, selectors, secondary_selector = gather_minerva_data(
            minerva_table, search_sets)
    x, y = search_sets[0][2].split(' ')
    for selector, found_match in zip(selectors, all_found_matches):
        data = {}
        found_all = [ found[0] for name, found in found_match.items() ]
        labels =  [ name for name, found in found_match.items() ] 
        print(found_all)
        # necessary because nested list
        for found, label in zip(found_all, labels):
            print(label)
            x_numbers = [ float(num[0]) for num in found ] 
            y_numbers = [ float(num[1]) for num in found ] 
            second_sel = [ label for i in range(len(x_numbers)) ]
            prim_sel = [ selector for i in range(len(x_numbers)) ]
            data[x] = x_numbers
            data[y] = y_numbers
            data[secondary_selector] = second_sel
            data["Selector"] = prim_sel
            print(data)
            #print(found_numbers) 
            #print(label)
            df_new = pd.DataFrame(data)
            df = df.append(df_new, ignore_index=True)
    print(df)
    print(secondary_selector)
    if secondary_selector:
        pp = sns.pairplot(data=df, x_vars=x, y_vars=y, hue=secondary_selector, kind='reg') 
    else:
        pp = sns.pairplot(data=df, x_vars=x, y_vars=y, kind='reg') 
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
    parser.add_argument("--pairplot", action='store_true', help="""Request that the 
    queries be summarised in a pairplot. This requires two magnitudes to be 
    requested.""")
    parser.add_argument("--stacked", action='store_true', help="""Request that the 
    queries be summarised in a stackedplot. Useful for comparing numbers of 
    taxa matching a query vs total.""")
    parser.add_argument("--boxplot", action='store_true', help="""Request that the 
    queries be summarised in a boxplot. Will summarise number the spread of 
    numbers of matching..........""")
    parser.add_argument("--normalize", action='store_true', help="""Attempts to 
    normalize data plotted.""")
    parser.add_argument("-f", "--file", required=False, help="""File in which 
    mutiple seperate queries and columns summaries are requested. Structured with 
    seperate queries with space, queries and columns seperated by tab, and 
    seperated searches by newline.""")

    args = parser.parse_args()

    if args.pairplot:
        with open(args.file) as queries:
            search_sets = [ search_set.strip().split('\t') for search_set in queries ]
        create_residual_summary(args.minervaTable, search_sets)
        sys.exit()

    if args.boxplot:
        with open(args.file) as queries:
            search_sets = [ search_set.strip().split('\t') for search_set in queries ]
        create_frequency_summary(args.minervaTable, search_sets)
        sys.exit()


    # Garbage temp code below
    if args.stacked or args.boxplot:
        if args.file:
            with open(args.file) as queries:
                found_series = defaultdict(list)
                found_series_set = []
                selector_set = []
                for each in queries:
                    found_set, selector, secondary_selection = _worker(args.minervaTable, 
                            each.strip().split('\t'))
                    print(secondary_selection)
                    print(selector)
                    if selector not in selector_set and selector_set:
                        cprint("appending", "red")
                        found_series_set.append(found_series)
                        found_series = defaultdict(list)
                        selector_set.append(selector)
                    elif not selector_set:
                        selector_set.append(selector)
                    # merge dicts
                    found_series = {**found_series, **found_set}
                found_series_set.append(found_series)
                print(len(found_series_set))
        if args.stacked:
            create_stacked_summary(found_series_set, secondary_selection)
        else:
            create_frequency_summary(found_series_set, selector_set, secondary_selection)

if __name__ == "__main__":
    main()
