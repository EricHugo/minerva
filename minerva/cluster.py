#!/usr/bin/env python

import sys
import subprocess
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

from minerva import diamondBlast

class clusterProteins():
    def __init__(self, result_tab, threads=1):
        self.blast_results = []
        df = pd.read_csv(result_tab, sep='\t')
        paths = df['Forward_path'].values.tolist()
        paths = paths + df['Reverse_path'].values.tolist()
        print(paths)
        self.blast_all(paths)
        self._format_blast()
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

    def _format_blast(self):
        results_noself = []
        all_genes = set()
        # remove self-self hits
        for blast_result in self.blast_results:
            if blast_result[0] == blast_result[1]:
                continue
            results_noself.append(blast_result)
            all_genes.add(blast_result[0])
            all_genes.add(blast_result[1])
        # fix order
        all_genes = list(all_genes)
        print(results_noself)
        print(all_genes)
        results_dict = defaultdict(list)
        results_dict['Gene'] = []
        for result in results_noself:
            results_dict['Gene'].append(result[0])
            for gene in all_genes:
                if result[1] == gene:
                    results_dict[gene].append(float(result[2]))
                else:
                    results_dict[gene].append(0)
        print(results_dict)
        df = pd.DataFrame.from_dict(results_dict)
        print(df.dtypes)
        print(df)
        return

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
        # blast cat file against db
        blast = diamondBlast(filename, filename)
        results = blast.perform_blast(filename,  '-f', '6', 'qseqid', 'sseqid',
                                      'evalue', 'pident', evalue="1e-2"
                                      ).strip().split('\n')
        for blast_result in results:
            self.blast_results.append(blast_result.split('\t'))
        print(self.blast_results)
        # return scores
        return self.blast_results

    def mcl_cluster(self):
        pass

    def assign_groups(self):
        pass
