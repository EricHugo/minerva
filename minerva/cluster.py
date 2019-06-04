#!/usr/bin/env python

import sys
import subprocess
import pandas as pd
from Bio import SeqIO

from minerva import diamondBlast

class clusterProteins():
    def __init__(self, result_tab, threads=1):
        df = pd.read_csv(result_tab, sep='\t')
        paths = df['Forward_path'].values.tolist()
        paths = paths + df['Reverse_path'].values.tolist()
        print(paths)
        self.blast_all(paths)
        return
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
        filename = 'all_neighbours.faa'
        handle = open(filename, 'w+')
        # cat all files to single file
        sequences = [SeqIO.parse(seq, 'fasta') for seq in files]
        for sequence in sequences:
            for seq in sequence:
                handle.write('>' + seq.id + '\n')
                handle.write(str(seq.seq) + '\n')
        handle.close()
        # makedb from cat file
        makedb = ['diamond', 'makedb', '--in', filename, '-d', filename]
        db = subprocess.run(makedb)
        if db.returncode > 0:
            print("Failed to makedb")
            return None
        blast = diamondBlast(filename, filename)
        blast_result = blast.perform_blast(filename, evalue="1e-2")
        print(blast_result)
        # blast cat file against db
        # return scores
        pass

    def mcl_cluster(self):
        pass

    def assign_groups(self):
        pass
