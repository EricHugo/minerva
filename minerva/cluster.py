#!/usr/bin/env python

import sys
import subprocess

from minvera import diamondBlast

class clusterProteins():
    def __init__(self, result_tab, threads=1):
        with open(result_tab) as f:
            self.result_tab = [entry.strip().split('\t') for entry in f]
        for n, entry in enumerate(self.result_tab[1]):
            if entry == 'Forward_path':
                fn = n
            elif entry == 'Reverse_path':
                rn = n
        self.forward_paths = {entry[fn - 4]: entry[fn] for entry in
                              self.result_tab}
        self.reverse_paths = {entry[rn - 4]: entry[rn] for entry in
                              self.result_tab}

    def blast_all(self, files):
        # cat all files to single file
        subprocess.run('cat'
        # makedb from cat file
        # blast cat file against db
        # return scores
        pass

    def mcl_cluster(self):
        pass

    def assign_groups(self):
        pass
