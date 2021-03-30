#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
from ftplib import FTP
from datetime import datetime
import logging
import sys
import mmap
import re
import os
import csv
import tarfile
import hashlib


class parseTaxonomy():
    def __init__(self, tax_path='taxonomy'):
        self.tax_path = tax_path
        self.download_taxonomy()
        self.tax_names = tax_path + '/names.dmp'
        self.tax_nodes = tax_path + '/nodes.dmp'

    def download_taxonomy(self):
        """Module downloads, unpacks, and verifies ncbi taxdump"""
        try:
            modified_time = datetime.fromtimestamp(os.path.getmtime(self.tax_path +
                '/nodes.dmp'))
            age = datetime.today() - modified_time
            if not age.days >= 2:
                return
        except FileNotFoundError:
            pass
        with FTP('ftp.ncbi.nih.gov') as ncbi:
            ncbi.login()
            ncbi.cwd('pub/taxonomy/')
            print("Downloading taxdump from NCBI...", file=sys.stderr, end='',
                    flush=True)
            ncbi.retrbinary('RETR taxdump.tar.gz', open('taxdump.tar.gz',
                'wb').write)
            ncbi.retrbinary('RETR taxdump.tar.gz.md5', open('taxdump.tar.gz.md5',
                'wb').write)
            print(" Done", file=sys.stderr)
            ncbi.quit()
        try:
            os.mkdir(self.tax_path)
        except FileExistsError:
            pass
        local_md5 = hashlib.md5(open('taxdump.tar.gz', 'rb').read()).hexdigest()
        with open('taxdump.tar.gz.md5') as md5:
            ncbi_md5 = md5.readline().split()[0]
            if not ncbi_md5 == local_md5:
                raise ValueError('Local md5 checksum does not equal value from ncbi')
        tar = tarfile.open('taxdump.tar.gz')
        tar.extractall(path=self.tax_path)
        os.remove('taxdump.tar.gz')
        os.remove('taxdump.tar.gz.md5')

    def find_taxid(self, name):
        """Given full scientific name: return taxid"""
        taxid = ''
        with open(self.tax_names, 'r') as names_file:
            for tax in names_file:
                if re.fullmatch(re.escape(name), tax.split('|')[1].strip()):
                    taxid = tax.split('|')[0].strip()
                    break
        return taxid

    def find_scientific_name(self, taxid):
        """Given taxid, return full scientific name"""
        with open(self.tax_names, 'r') as names_file:
            for tax in names_file:
                if re.fullmatch(taxid, tax.split('|')[0].strip()) and \
                   re.fullmatch('scientific name', tax.split('|')[3].strip()):
                    name = tax.split('|')[1].strip()
                    break
        return name

    def parse_taxa(self, taxid):
        """Will locate parent of given taxid, recursively loop to find full i
        taxonomic list, returns list of tuple(rank, taxid) pairs"""
        taxids = taxid.split()
        ranks = []
        with open(self.tax_nodes, 'r+b') as node_file:
            m = mmap.mmap(node_file.fileno(), 0, prot=mmap.PROT_READ)
            while True:
                # dict comp to quickly find match first field and retrieve second
                parent = [ (taxids.append(str(node, "utf-8").split('|')[1].strip()),
                        ranks.append(str(node, "utf-8").split('|')[2].strip()))
                        for node in iter(m.readline, bytes())
                        if re.fullmatch(str(taxids[-1]), str(node, "utf-8").
                            split('|')[0].strip()) ]
                m.seek(0)
                if taxids[-1] == '1':
                    taxids.pop()
                    ranks[0] = "Name"
                    break
        return list(zip(ranks, taxids))

