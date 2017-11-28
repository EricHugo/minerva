#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from micomplete import calcCompleteness
from micomplete import parseSeqStats
#from micomplete import micomplete
from contextlib import contextmanager
from distutils import spawn
from collections import defaultdict
from ftplib import FTP
from datetime import datetime
from termcolor import cprint
import sys
import mmap
import re
import os
import csv
import argparse
import shutil
import multiprocessing as mp
import subprocess
import tarfile
import hashlib
import subprocess

# import dev version of miComplete
import importlib.util
spec = importlib.util.spec_from_file_location("micomplete", 
"/home/hugoson/git/micomplete/micomplete/micomplete.py")
micomplete = importlib.util.module_from_spec(spec)
spec.loader.exec_module(micomplete)

# import dev minerva modules
sys.path.append('/home/hugoson/git/minerva')
try:
    from minerva import parseTaxonomy, parseMapFile
except ImportError:
    from parse_taxonomy import parseTaxonomy
    from parse_minerva_map import parseMapFile


def _worker(fasta, seqType, name, hmm, q, gen_directory, evalue=1e-20, 
            crispr=False, outfile=None):
    tax = parseTaxonomy()
    if seqType == "faa":
        faa = fasta
        return #not supported for now
    elif seqType == "fna":
        faa = micomplete.create_proteome(fasta)
        return #not supported for now
    elif re.match("(gb.?.?)|genbank", seqType):
        name = get_gbk_feature(fasta, 'organism')
        print("Working on %s" % name)
        faa = micomplete.extract_gbk_trans(fasta, re.sub('\)|\(|\{|\}|\[|\]|\/|\/', 
            '', name) + '.faa')
        fna = get_contigs_gbk(fasta, re.sub('\/', '', name))
        # if there is no translation to extract, get contigs and use prodigal
        # find ORFs instead
        if os.stat(faa).st_size == 0:
            os.remove(faa)
            faa = micomplete.create_proteome(fna, re.sub('\/', '', name))
        # find CRISPRs
        if crispr:
            c_out, crispr_bool = find_CRISPRs(fna, re.sub('\/', '', name))
        else:
            c_out, crispr_bool = "-", "-"
        # grab size and GC stats
        pseqs = parseSeqStats(fasta, name, seqType)
        seqLength, alllengths, GC = pseqs.get_length()
        taxid = tax.find_taxid(name)
        if taxid:
            lineage = tax.parse_taxa(taxid)
            taxonomy = { rank: tax.find_scientific_name(taxid) for rank, taxid 
                    in lineage }
        else:
            taxonomy = {}
    else:
        raise TypeError('Sequence type needs to be specified as one of faa/fna/gbk')
    baseName = os.path.basename(fasta).split('.')[0]
    if not name:
        name = baseName
    gene_matches = get_matches(faa, name, hmm, evalue)
    for gene, match in gene_matches.items():
        gene_matches[gene].append(extract_protein(faa, gene))
    # make an entry of empty results  
    if not gene_matches:
        gene_matches['-'].append(['-', '-'])
    compile_results(name, gene_matches, taxid, taxonomy, fasta, seqType, faa, q, 
            seqLength, GC, gen_directory=gen_directory, crispr_bool=crispr_bool,
            c_out=c_out)
    return 

def compile_results(name, gene_matches, taxid, taxonomy, fasta, seqType, faa, q,
        seqLength, GC, gen_directory="protein_matches", crispr_bool="-", c_out="-"):
    for gene, match in gene_matches.items():
        #print(match)
        result = {}
        result['name'] = name
        result['match'] = match[0][0]
        result['evalue'] = match[0][1]
        result['gene'] = gene
        result['taxid'] = taxid
        result.update(taxonomy)
        # try to output the gene-file, pass if no genes were found
        try:
            result['gene_path'] = os.path.abspath(gen_directory + '/' + 
                    match[1][0].name + ".faa")
            q.put(match[1])
        except IndexError:
            pass
        result['genome_path'] = os.path.abspath(fasta)
        result['proteome_path'] = os.path.abspath(faa)
        result['crispr'] = str(crispr_bool)
        result['crispr_path'] = c_out
        result['genome_length'] = str(seqLength)
        result['gc-content'] = str(GC)
        # put result dict in queue for listener
        q.put(result)
        # also put the seq-object
    return

def _listener(q, headers, outfile='-', gen_directory="protein_matches"):
    """Process to write results in a thread safe manner"""
    try:
        os.mkdir("protein_matches")
    except FileExistsError:
        pass
    with open_stdout(outfile) as handle:
        while True:
            out_object = q.get()
            if out_object == "done":
                break
            elif type(out_object) is dict:
                for header in headers:
                    try:
                        handle.write(out_object[header.lower()] + '\t')
                    except KeyError:
                        handle.write('-' + '\t')
                handle.write('\n')
            if type(out_object) is list:
                # write aminoacid sequence
                for seq_object in out_object:
                    try:
                        SeqIO.write(seq_object, gen_directory + '/' +
                                seq_object.name + ".faa", "fasta")
                    except (IOError, AttributeError):
                        cprint("Unable to output sequence: " + seq_object, "red", 
                                file=sys.stderr)
            elif type(out_object) is tuple:
                for head in out_object:
                    handle.write(head + "\t")
                handle.write("\n")
            else:
                continue
    return

@contextmanager
def open_stdout(outfile):
    """Function dynamically allows printing to file or stdout when used as open 
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

def get_matches(faa, name, hmm, evalue=1e-20):
    """Retrieves markers in hmm from the given proteome"""
    comp = calcCompleteness(faa, re.sub('\/', '', name), hmm, evalue)
    foundmatches, dupmatches, totl = comp.get_completeness()
    print(name + ': ' + str(foundmatches) )
    # ensure unique, best match for each hmm
    # but allow infinite gene matches for each hmm 
    gene_matches = defaultdict(list)
    try:
        for matche, match in foundmatches.items():
            for gene, evalue in match:
                if gene not in gene_matches:
                    gene_matches[gene].append([matche, evalue])
                elif float(evalue) < float(gene_matches[gene][0][1]):
                    gene_matches[gene].pop()
                    gene_matches[gene].append([matche, evalue])
    except AttributeError:
        print(name + " threw AttributeError")
    return gene_matches

def extract_protein(faa, gene_tag):
    """Extract the protein sequences from a given proteome given a dict of 
    gene names"""
    #rayt_gene_list = [ gene for gene, rayt in gene_RAYT.items() ]
    protein_list = [ seq for seq in SeqIO.parse(faa, "fasta") if seq.id in
            gene_tag ]
    #print(protein_list)
    return protein_list

def get_gbk_feature(handle, feature_type):
    """Get specified organism feature from gbk file"""
    input_handle = open(handle, mode='r')
    for record in SeqIO.parse(input_handle, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                value = ''.join(feature.qualifiers[feature_type])
                break
    return value

def get_contigs_gbk(gbk, name):
    """Extracts all sequences from gbk file, returns filename"""
    handle = open(gbk, mode='r')
    if not name:
        name = os.basename(gbk).split('.')[0]
    out_handle = open(name, mode='w')
    for seq in SeqIO.parse(handle, "genbank"):
        out_handle.write(">" + seq.id + "\n")
        out_handle.write(str(seq.seq) + "\n")
    out_handle.close()
    return name

def init_results_table(q, outfile=None):
    headers = [
            'Match',
            'Gene',
            'Evalue',
            'Name',
            'TaxID',
            'Superkingdom',
            'Phylum',
            'Class',
            'Order',
            'Family',
            'Genus',
            'Species',
            'Gene_path',
            'Genome_path',
            'Proteome_path',
            'CRISPR',
            'CRISPR_path',
            'Genome_length',
            'GC-content'
            ]
    if outfile and not outfile == '-':
        headers = tuple(headers)
        q.put(headers)
    else:
        with open_stdout(outfile) as handle:    
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(headers)
    return headers

def find_CRISPRs(fna, name=None):
    """Uses CRISPR recognition tool (CRT) to find putative CRISPR repeats in 
    given fna file, returns out-table, as well as True / False for positive / 
    negative result."""
    input_handle = open(fna, mode='r')
    if name:
        output_handle = name + ".cout"
    else:
        baseName = os.path.basename(gbkfile).split('.')[0]
        outfile = baseName + ".cout"
        output_handle = outfile
    run_fmt = [ 'java', '-cp', 'CRT1.2-CLI.jar', 'crt', fna, output_handle ]
    return_proc = subprocess.run(run_fmt, stdout=subprocess.DEVNULL)
    errcode = return_proc.returncode
    if errcode > 0:
        print(return_proc)
        raise RuntimeError('Something went wrong with CRISPR finder for file %s'
                % fna)
    with open(output_handle) as out:
        for line in out:
            if re.match("CRISPR [0-9]", line):
                result = True
                break
            if re.match("No CRISPR elements were found", line):
                result = False
                break
    return output_handle, result

def main():
    parser = argparse.ArgumentParser(description="""For a given set of genomes
            investigate presence of RAYT via HMMs. Attempt to categorize mathces
            and investigate flanking genes and 16mers.""")

    parser.add_argument("fastaList", help="""Sequence(s) along with type (fna, 
            faa, gbk) provided in a tabular format""")
    parser.add_argument("-o", "--outfile", required=False, default="-",
            help="""Filename to save results. Otherwise prints to stdout.""")
    parser.add_argument("--gendir", required=False, default="protein_matches",
            help="""Directory in which to store matched protein sequences""")
    parser.add_argument("--hmms", required=False, default=False,
            help="""Specifies a set of HMMs to be used for completeness check 
                        or linkage analysis""")
    parser.add_argument("--taxa", required=False, help="""Query specific taxonomic
            group, requires a csv of the appropriate group from the NCBI genome
            browser""")
    parser.add_argument("--glist", required=False, help="""Genome list in csv 
            format from the NCBI genome browser, required with '--taxa' 
                        argument""")
    parser.add_argument("--crispr", required=False, action='store_true',
            help="""Flag to attempt to assign CRISPR systems within examined 
            genomes using CRISPR Recognition Tool (CRT).""")
    parser.add_argument("--summary", required=False, help="""Attempts to provide 
            a summary of pre-exist results file. Provide a file of column(s) to 
            be summarised and optionally a selection column with a string to 
            be matched within the selection column""")
    parser.add_argument("--threads", required=False, default=1, type=int,
            help="""Number of threads to be used concurrently""")
    args = parser.parse_args()

    try:
        assert shutil.which('hmmsearch')
    except AssertionError:
        raise RuntimeError('Unable to find hmmsearch in path')
    except AttributeError:
        try:
            assert spawn.find_executable('hmmsearch')
        except AssertionError:
            raise RuntimeError('Unable to find hmmsearch in path')

    # Initialise taxdump, threadsafety
    parseTaxonomy()

    with open(args.fastaList) as seq_file:
        inputSeqs = [ seq.strip().split('\t') for seq in seq_file ]
    manager = mp.Manager()
    q = manager.Queue()
    headers = init_results_table(q, args.outfile)
    pool = mp.Pool(processes=int(args.threads) + 1)
    # init listener here
    listener = pool.apply_async(_listener, (q, headers, args.outfile, args.gendir))
    
    jobs = []
    for i in inputSeqs: 
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(_worker, (i[0], i[1], i[2], args.hmms, q, 
            args.gendir), {"crispr":args.crispr})
        jobs.append(job)
    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    listener.get()
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
