============
**minerva**
============


- Eric Hugoson (eric@hugoson.org / hugoson@evolbio.mpg.de / @e_hugoson)


Introduction
--------------
With the increasing availability of genomic data in online databases such as NCBI
there is a need to be able to systematically explore the genomic diversity.
Particularly, a large subset of genes' functions are partly or entirely unknown.
Large sets of genomic data allows exploration of gene function on a different level
from before, allowing inference of function from differences between which
genomes by presence/absence.

minvera enables the systematic exploration of microbial genomes for
presence/absence of a single- or a set of gene(s) given as a Hidden Markov Model.
While identifying presence/absence of genes minerva also identifies other
relevant statistic of each genome, permitting underlying patterns to be identified and
ultimately for potential functions to be inferred.

minerva is still in an early stage of development. As such, bugs are very much
to be expected.

Description
--------------
minvera is a software designed to build a database of parsed genomes along with
their relevant attributes (e.g. taxonomy) with presence/absence of a gene or
genes of interest.

Dependencies
--------------
Python (>=3.5)


External software
^^^^^^^^^^^^^^^^^^^
Executables should be available in the user's ``$PATH``.

HMMER3
"""""""""""""""""
HMMER: biosequence analysis using profile hidden Markov models, by Sean Eddy and coworkers. Tested with v. 3.3. Available from <http://hmmer.org/>.

prodigal
""""""""""""""""
A gene prediction software by Doug Hyatt. Tested with v. 2.6.3. Download at:
<https://github.com/hyattpd/Prodigal>

Python libraries
^^^^^^^^^^^^^^^^^^^
If built from the package these will be installed automatically, otherwise can easily be installed seperately using ``pip`` (`Install pip <https://pip.pypa.io/en/stable/installing/>`_).

Required
""""""""""""""""""

- Biopython (>= 1.70) (``$ pip install biopython``)
- Numpy (>= 1.13.1) (``$ pip install numpy``)
- Matplotlib (>= 2.0.2) (``$ pip install matplotlib``)
- miComplete (>= 1.0.0) (``$ pip install miComplete``)

Installation
--------------

Python3
^^^^^^^

``$ pip install PLACEHOLDER``

Usage
--------------

Download Genbank Files
^^^^^^^^^^^^^^^^^^^^^^^
How do download genbank files from ncbi, specific groups. With example(s). Also link a more comprehensive guide.

Positional arguments
^^^^^^^^^^^^^^^^^^^^^^^
::

  genomes               Sequence(s) along with type (fna, faa, gbk) provided
                        in a tabular format
  hmms                  File containing one or more HMM of genes to be
                        identified in given genomes.



The ``genomes`` file has to contain a path per line (relative or absolute) to a genbank file: ::

   /seq/genomic_sources/legionella_pneumophila.gbff
   /seq/genomic_sources/coxiella_burnetii.gbff
   /seq/genomic_sources/e_coli.gbff
   (...)

The ``hmms`` has to be a file containing one or more Hidden Markov Models (HMM) produced by HMMER. These can be produced via HMMER directly from aligned sequence files or obtained from many online resources such as eggNOG.

Optional arguments
^^^^^^^^^^^^^^^^^^^^^^^^

  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Filename to save results. Otherwise prints to stdout.
  --gendir GENDIR       Specify directory in which to store matched protein
                        sequences.
  --datadir DATADIR     Specify directory in which to store produced genomes,
                        proteomes, and HMMER results.
  --taxa TAXA           Query specific taxonomic group, requires a csv of the
                        appropriate group from the NCBI genome browser
  --rank RANK           Rank (e.g. "order, family, genus) of taxa specified
                        in --taxa argument. Required with the --taxa argument
  --crispr              Flag to attempt to assign CRISPR systems within
                        examined genomes using CRISPR Recognition Tool (CRT).
  --no_neighbours       Flag to skip identifying nearest neighbour up- and
                        downstream of matched gene.
  --threads THREADS     Number of threads to be used concurrently.
  --db DB               Diamond database to blast neighbour proteins against
                        for product description if labeled hypothetical in
                        gbk.



Examples
^^^^^^^^^^^^^^^^^^^^^^^^
In a folder containing one or several FASTA files with '.gbk' extensions, create a list of genbank files. Here it is best to avoid relative paths unless you know you will be running minerva from the same relative directory. A correctly formatted input file can be created via a simple ``find`` command::

   $ find $(realpath .) -maxdepth 1 -type f -name "*.gbff" > minerva_in.list

Note that the files names must contain a valid suffix. All variations of genbank suffixes (e.g. .gb, .gbk, gbff) should be valid.

Sequence list file, ``minerva_in.list``::

   $ cat minerva_in.list
    /seq/genomic_sources/legionella_longbeachae.gbk
    /seq/genomic_sources/coxiella_burnetii.gbk
    /seq/genomic_sources/coxiella-like_endosymbiont.gbk
    /seq/genomic_sources/e_coli.gbk

Additionally a HMM file for one or many genes is required. This one was aquired from EMBL's orthology database, eggNOG, in the Legionellales dataset ``http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/118969/118969_hmms.tar.gz``::

    1JER9.hmm

Example 1 - Basic presence/absence
""""""""""""""""""""""""""""""""""
Here we are only looking to go through all genomes in the created ``minerva_in.list`` to see where the gene !!! is present and where it is not. Therefore we also add the ``--noneighbours`` argument so minerva will not attempt to identify neighbour genes of matched genes::

    $ minerva minerva_in.list 1JER9.hmm --noneighbours

This produces results in a tabular format directly to stdout which can be difficult to read. It is therefore better to redirect it to a file and viewed in a spreadsheet reader::

    $ minerva minerva_in.list 1JER9.hmm --noneighbours > minerva_results.tsv
    $ libreoffice --calc minerva_results.tsv


Example 2 - Presence/absence with neighbourhood
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
If we wish to know more about the genes, minerva is capable of locating and describing the neighbourhood of each matched gene by leaving out the ``--noneighbours`` flag. We can also organise the output directory by telling minerva where to store produced files::

    $ minerva minerva_in.list 1JER9.hmm -o minerva_results.tsv --gendir matched_genes --datadir produced_genomes

This will leave our working directory empty save for the results file ``minerva_results.tsv``, which now also contains information regarding the neighbourhood of matched genes. However, the runtime increased significantly due to needing to map the neighbourhoods, therefore we may wish to run multithreaded::

    $ minerva minerva_in.list 1JER9.hmm -o minerva_results.tsv --gendir matched_genes --datadir produced_genomes --threads 4

