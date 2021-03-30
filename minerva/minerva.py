#!/usr/bin/env python

from __future__ import print_function
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
import logging
import logging.handlers
import subprocess
from Bio import SeqIO
from micomplete import calcCompleteness, parseSeqStats, micomplete
from contextlib import contextmanager
from distutils import spawn
from collections import defaultdict
from itertools import zip_longest
from ftplib import FTP
from datetime import datetime
from termcolor import cprint
from pathlib import Path


# import dev minerva modules
sys.path.append('/home/hugoson/git/minerva')
try:
    from minerva import parseTaxonomy, parseMapFile, findGeneNeighbourhood, diamondBlast, clusterProteins
except ImportError:
    from parse_taxonomy import parseTaxonomy
    from parse_minerva_map import parseMapFile
    from find_gene_neighbourhood import findGeneNeighbourhood
    from diamondblast import diamondBlast
    from cluster import clusterProteins

ILLEGAL_CHARACTERS = "|><()[]{}=*/"
NEIGHBOUR_DIR = "neighbour_proteins"
try:
    os.mkdir(NEIGHBOUR_DIR)
except FileExistsError:
    pass

def _worker(genome, seqType, raw_name, hmm, q, gen_directory, datadir, tax, evalue=1e-30,
            crispr=False, blast_db=None, taxa=None, neighbourhood=True):
    #tax = parseTaxonomy()
    if seqType == "faa":
        faa = genome
        return #not supported for now
    elif seqType == "fna":
        faa = micomplete.create_proteome(genome)
        return #not supported for now
    elif re.match("(gb.?.?)|genbank", seqType):
        for feature in get_gbk_feature(genome, 'features'):
            if feature.type == "source" and not raw_name:
                raw_name = ''.join(feature.qualifiers['organism'])
                break
        taxid = tax.find_taxid(raw_name)
        if taxid:
            lineage = tax.parse_taxa(taxid)
            taxonomy = {rank: tax.find_scientific_name(taxid) for rank, taxid
                        in lineage}
        else:
            taxonomy = {}
        print(taxonomy)
        name = raw_name
        # check if desired taxa, else skip
        if taxa[0]:
            if not taxa[0].lower() == taxonomy[taxa[1].lower()].lower():
                return
        for i in ILLEGAL_CHARACTERS:
            name = re.sub(re.escape(i), '', name)
        name = re.sub(' ', '_', name)

        # get strain name to prevent collisions
        for feature in get_gbk_feature(genome, 'features'):
            if feature.type == "source":
                try:
                    strain = ''.join(feature.qualifiers['strain'])
                    try:
                        strain = ''.join(strain.split(':')[1])
                    except IndexError:
                        pass
                    break
                except KeyError:
                    strain = ''.join(feature.qualifiers['db_xref'])
        name = raw_name + '_' + strain
        for i in ILLEGAL_CHARACTERS:
            name = re.sub(re.escape(i), '', name)
        name = re.sub(' ', '_', name)

        print("Working on %s" % raw_name)
        # remove illegal characters
        # also swap spaces for underscores
        faa_name = os.path.join(datadir, name)
        faa = micomplete.extract_gbk_trans(genome, faa_name + '.faa')
        fna = get_contigs_gbk(genome, datadir, re.sub('\/', '', name + '.fna'))
        # if there is no translation to extract, get contigs and use prodigal
        # find ORFs instead
        if not os.path.isfile(faa) or os.stat(faa).st_size == 0:
            print("prodigal")
            cprint(faa_name, "red")
            try:
                os.remove(faa)
            except FileNotFoundError:
                pass
            faa = micomplete.create_proteome(fna, faa_name + '.faa')
        # find CRISPRs
        if crispr:
            c_out, c_bool = find_CRISPRs(fna, re.sub('\/', '', name))
        else:
            c_out, c_bool = "-", "-"
        # grab size and gc_content stats
        pseqs = parseSeqStats(genome, name, seqType)
        seq_length, _, gc_content = pseqs.get_length()
    else:
        raise TypeError('Sequence type needs to be specified as one of '
                        'faa/fna/gbk')
    basename = os.path.basename(genome).split('.')[0]
    if not name:
        name = basename
    gene_matches = get_matches(faa, faa_name, hmm, evalue)
    print(gene_matches)
    print(neighbourhood)
    if gene_matches and neighbourhood:
        neighbour_find = findGeneNeighbourhood(genome, faa, name, seq_length,
                                               gene_matches)
        neighbours = neighbour_find.find_minimum_distance()

        # get product names if possible
        print(neighbours)
        if neighbours:
            for gene, neighbourhood in neighbours.items():
                print(gene)
                neighbours = get_neighbour_products(faa, genome, neighbourhood, gene, q, blast_db)
                gene_matches[gene].append(extract_protein(faa, gene))
                gene_matches[gene].append(neighbours)
                print(gene_matches)
    elif gene_matches:
        pass
    else:
        gene_matches['-'].append(['-', '-'])
    # make an entry of empty results
    compile_results(name, gene_matches, taxid, taxonomy, genome, seqType, faa, q,
                    seq_length, gc_content, gen_directory=gen_directory, c_bool=c_bool,
                    c_out=c_out)
    return

def get_neighbour_products(faa, fasta, neighbours, gene, q, blast_db):
    # neighbour_OGs = []
    for each in ['forward', 'reverse']:
        # if no neighbour for this particular direction, then append no product
        if not neighbours[each + '_neighbour']:
            continue
        neighbour_product = None
        for prod in find_gbk_product(fasta, target=neighbours[each + '_neighbour'],
                                     unique=True):
            cprint(prod, "red")
            neighbour_faa = extract_protein(faa, neighbours[each + '_neighbour'],
                                            write=True)
            if ''.join(prod) == "hypothetical protein" and blast_db:
                continue
            neighbour_product = prod
            print(neighbour_faa)
        if not neighbour_product and blast_db:
            cprint(gene, "magenta")
            print(neighbours[each + '_neighbour'])
            print(faa)
            neighbour_faa = extract_protein(faa, neighbours[each + '_neighbour'],
                                            write=True)
            # new func that takes extracted protein and diamond
            # against uniprot => extracts protein name
            blast = diamondBlast(blast_db, '_')
            cprint(neighbour_faa, "magenta")
            try:
                blast.perform_blast(neighbour_faa, '--outfmt', '6',
                                    'sseqid', 'evalue', 'bitscore', 'ppos',
                                    'salltitles', evalue='1e-10')
            except (TypeError, IndexError):
                out = [fasta, faa, neighbours[each + '_neighbour']]
                cprint(out, "red")
            try:
                neighbours[each + '_product'] = ''.join(blast.get_protein_name())
            except TypeError:
                cprint(blast.get_protein_name(), "red")
                neighbours[each + '_product'] = "hypothetical protein"
        else:
            neighbours[each + '_product'] = "-"
        neighbours[each + '_path'] = os.path.abspath(neighbour_faa)
        #neighbour_OG = get_matches(neighbour_faa, "test", '/seq/databases/eggNOG/bactNOG/all_bactNOG.hmm')
        #neighbours[each + '_OG'] = list(neighbour_OG.values())[0][0][0]
        neighbours[each + '_OG'] = "NA"
    print(neighbours)
    return neighbours

def compile_results(name, gene_matches, taxid, taxonomy, fasta, seqType, faa, q,
        seq_length, gc_content, gen_directory="protein_matches", c_bool="-", c_out="-"):
    for gene, match in gene_matches.items():
        print(match)
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
        result['crispr'] = str(c_bool)
        result['crispr_path'] = c_out
        result['genome_length'] = str(seq_length)
        result['gc-content'] = str(gc_content)
        try:
            result['forward_neighbour'] = match[2]['forward_neighbour']
            result.update(match[2])
            #result['reverse_neighbour'] = match[2][1][0]
            #if match[2][0][1]:
            #    result['forward_distance'] = str(match[2][0][1])
            #if match[2][1][1]:
            #    result['reverse_distance'] = str(match[2][1][1])
            #result['forward_product'] = match[2][2]
            #result['reverse_product'] = match[2][3]
            #print(match[2][3])
        except IndexError:
            pass
        # put result dict in queue for listener
        q.put(result)
    return


def _configure_logger(q, name, level=logging.WARNING):
    qh = CustomQueueHandler(q)
    logger = logging.getLogger(name)
    logformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    qh.setFormatter(logformatter)
    logger.addHandler(qh)
    logger.setLevel(level)
    return logger


def _listener(q, headers, outfile='-', gen_directory="protein_matches"):
    """Process to write results in a thread safe manner"""
    try:
        os.mkdir(gen_directory)
    except FileExistsError:
        pass
    with open_stdout(outfile) as handle:
        while True:
            out_object = q.get()
            if out_object == "done":
                break
            elif isinstance(out_object, dict):
                for header in headers:
                    try:
                        handle.write(out_object[header.lower()] + '\t')
                    except (KeyError, TypeError):
                        handle.write('-' + '\t')
                handle.write('\n')
            if isinstance(out_object, list):
                # write aminoacid sequence
                for seq_object in out_object:
                    try:
                        SeqIO.write(seq_object, gen_directory + '/' +
                                    seq_object.name + ".faa", "fasta")
                    except (IOError, AttributeError):
                        cprint("Unable to output sequence: " + seq_object, "red",
                               file=sys.stderr)
            elif isinstance(out_object, tuple):
                for head in out_object:
                    handle.write(head + "\t")
                handle.write("\n")
            else:
                continue
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

def get_matches(faa, name, hmm, evalue=1e-30):
    """Retrieves markers in hmm from the given proteome"""
    comp = calcCompleteness(faa, name, hmm, evalue=evalue, bias=0.3, best_domain=1e-5)
    foundmatches, _, _ = comp.get_completeness()
    print(name + ': ' + str(foundmatches))
    # ensure unique, best match for each hmm
    # but allow infinite gene matches for each hmm
    gene_matches = defaultdict(list)
    try:
        for matche, match in foundmatches.items():
            for gene, _, _, _, evalue in match:
                if gene not in gene_matches:
                    gene_matches[gene].append([matche, evalue])
                elif float(evalue) < float(gene_matches[gene][0][1]):
                    gene_matches[gene].pop()
                    gene_matches[gene].append([matche, evalue])
    except AttributeError:
        print(name + " threw AttributeError")
    return gene_matches

def extract_protein(faa, gene_tag, write=None):
    """Extract the protein sequences from a given proteome given a dict of
    gene names"""
    protein_list = [seq for seq in SeqIO.parse(faa, "fasta") if seq.id in
                    gene_tag]
    if write:
        print(protein_list)
        for protein in protein_list:
            SeqIO.write(protein, NEIGHBOUR_DIR + '/' + protein.name + '.faa',
                        'fasta')
            return NEIGHBOUR_DIR + '/' + protein.name + '.faa'
    return protein_list

def get_gbk_feature(handle, feature_type):
    """Get specified organism feature from gbk file"""
    input_handle = open(handle, mode='r')
    value = None
    for record in SeqIO.parse(input_handle, "genbank"):
        for feature in getattr(record, feature_type):
            #print(feature)
            yield feature
    return value

def get_contigs_gbk(gbk, dest=None, name=None):
    """Extracts all sequences from gbk file, returns filename"""
    handle = open(gbk, mode='r')
    if not name:
        name = os.path.basename(gbk).split('.')[0]
    if dest:
        try:
            name = os.path.join(dest, name)
        except FileNotFoundError:
            os.mkdir(dest)
    out_handle = open(name, mode='w')
    for seq in SeqIO.parse(handle, "genbank"):
        out_handle.write(">" + seq.id + "\n")
        out_handle.write(str(seq.seq) + "\n")
    out_handle.close()
    return name

def find_gbk_product(gbk, target=True, unique=False):
    """Matches target locus tag in given .gbk file. Yields any matched
    product. Yields all products in gbk if no target given"""
    # hacky index
    n = 0
    if unique and n > 0:
        return
    try:
        for feature in get_gbk_feature(gbk, "features"):
            if feature.type == 'CDS':
                if ''.join(feature.qualifiers['locus_tag']) == target:
                    n += 1
                    yield feature.qualifiers['product']
    except KeyError:
        return None
    return None

def init_results_table(q, outfile=None):
    headers = (
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
        'GC-content',
        'Forward_neighbour',
        'Forward_distance',
        'Forward_product',
        'Forward_OG',
        'Forward_path',
        'Reverse_neighbour',
        'Reverse_distance',
        'Reverse_product',
        'Reverse_OG',
        'Reverse_path',
        )
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
    run_fmt = ['java', '-cp', 'CRT1.2-CLI.jar', 'crt', fna, output_handle]
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
    parser = argparse.ArgumentParser(
        description="""Identification of genes within a given set of microbial 
                genomes. minerva maps presence/absence of gene families, and 
                assists in identifying patterns underlying the 
                presence/absence.""",
        epilog="""Report issues and bugs to the issue tracker at
                https://github.com/EricHugo/minerva or directly to 
                eric@hugoson.org""",
        prog='minerva')

    parser.add_argument("genomes", help="""Sequence(s) along with type (fna,
            faa, gbk) provided in a tabular format""")
    parser.add_argument("hmms", help="""File containing one or more HMM of
            genes to be identified in given genomes.""")
    parser.add_argument("-o", "--outfile", required=False, type=str, default="-",
            help="""Filename to save results. Otherwise prints to stdout.""")
    parser.add_argument("--gendir", required=False, default="protein_matches",
            help="""Specify directory in which to store matched protein sequences.""")
    parser.add_argument("--datadir", required=False, default="data",
            help="""Specify directory in which to store produced genomes, proteomes,
            and HMMER results.""")
    parser.add_argument("--taxa", required=False, help="""Query specific taxonomic
            group, requires a csv of the appropriate group from the NCBI genome
            browser""")
    parser.add_argument("--rank", required=False, help="""Rank (e.g. "order, 
            family, genus) of taxa specified in --taxa argument.
            Required with the --taxa argument""")
    parser.add_argument("--crispr", required=False, action='store_true',
            help="""Flag to attempt to assign CRISPR systems within examined
            genomes using CRISPR Recognition Tool (CRT).""")
    parser.add_argument("--noneighbours", required=False, action='store_false',
            help="""Flag to skip identifying nearest neighbour up- and downstream
            of matched gene.""")
    parser.add_argument("--clusterneighbours", required=False, action='store_true',
            help="""Flag to cluster neighbours once identified. This defines
            clusters of orthologous genes to help identify related neighbourhoods.""")
    parser.add_argument("--threads", required=False, default=1, type=int,
                        help="""Number of threads to be used concurrently.""")
    parser.add_argument("--db", required=False, default=None, help="""Diamond
                        database to blast neighbour proteins against for
                        product description if labeled hypothetical in gbk.""")
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

    if args.taxa and not args.rank:
        raise RuntimeError('Taxa argument also requires specifying rank argument')

    # Initialise taxdump, threadsafety
    tax = parseTaxonomy()

    with open(args.genomes) as seq_file:
        input_seqs = [seq.strip().split('\t') for seq in seq_file if not
                      re.match('#|\n', seq)]

    manager = mp.Manager()
    q = manager.Queue()
    headers = init_results_table(q, args.outfile)
    pool = mp.Pool(processes=int(args.threads) + 1)
    # init listener here
    listener = pool.apply_async(_listener, (q, headers, args.outfile, args.gendir))
    jobs = []
    if args.datadir:
        os.makedirs(args.datadir, exist_ok=True)
    for i in input_seqs:
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(_worker, (i[0], i[1], i[2], args.hmms, q,
                                         args.gendir, args.datadir, tax),
                                         {"crispr":args.crispr,
                                          "blast_db":args.db,
                                          "taxa": (args.taxa, args.rank),
                                          "neighbourhood": args.noneighbours})
        jobs.append(job)
    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    listener.get()
    pool.close()
    pool.join()
    # finally cluster neighbours
    if args.noneighbours and args.clusterneighbours:
        clustering = clusterProteins(args.outfile, threads=args.threads)
        out = clustering.mcl_cluster()
        clusters = clustering.assign_groups(out)
        clustering.attribute_COGs(clusters)
    return

if __name__ == "__main__":
    main()
