# PhaMMseqs

The PhaMMseqs package facilitates pham assembly using [MMseqs2](https://www.mmseqs.com). Default parameters have
been carefully tuned for rapid, accurate exploration of the bacteriophage protein sequence space.

# Conda installation

The easiest way to install the phammseqs package and its dependencies is through the Anaconda/Miniconda package manager: 

    conda create -n phammseqs-env python=3.9 -y && conda activate phammseqs-env
    conda install -c bioconda -c conda-forge mmseqs2=13.45111 clustalo -y
    pip3 install phammseqs 
    
# Conda installation (Apple Silicon)

Macintosh computers purchased in the last couple of years no longer have Intel CPUs, instead sporting some flavor of Apple 
Silicon. These processors have a different architecture (arm64) that is not natively compatible with most of the recipes 
in the bioconda or conda-forge channels. We can run these programs using Apple's Rosetta2 emulator by modifying the conda
installation [as indicated here](https://github.com/Haydnspass/miniforge#rosetta-on-mac-with-apple-silicon-hardware).

    CONDA_SUBDIR=osx-64 conda create -n phammseqs-env python=3.9 -y && conda activate phammseqs-env
    conda env config vars set CONDA_SUBDIR=osx-64
    conda install -c bioconda -c conda-forge mmseqs2=13.45111 clustalo -y
    pip3 install phammseqs
    
After the second command you may be prompted to reactivate the environment for changes to take effect. This is easily achieved
by running this sequence of commands before you install any packages into the environment:

    conda deactivate && conda activate phammseqs-env

# Manual installation

If you don't have some flavor of conda available (and don't want to install it...) you may follow the instructions
[here](https://github.com/soedinglab/mmseqs2#installation) to manually install `mmseqs`. An optional dependency,
`clustalo` can be manually installed following the instructions [here](http://www.clustal.org/omega/). 
Most modern operating systems also ship with Python3, the programming language used to develop this package, and 
required to run it. However, if your system does not have Python 3.7 or higher, you will need to obtain it 
[here](https://www.python.org/downloads/).

Once all that is done, you can obtain the phammseqs package from PyPI using pip:

    pip3 install phammseqs

# Using PhaMMseqs at the command line

If you installed phammseqs and its dependencies using either of the Conda approaches, you will need to activate the
environment before using `phammseqs` (substitute `phammseqs-env` with whatever you named the environment):

    conda activate phammseqs-env

You can invoke `phammseqs` with the `-h` option to print the help menu:

    phammseqs -h

Which should print something like:

    usage: phammseqs [-h] [--identity] [--coverage] [--evalue] [--sensitivity] [--cluster-mode] [--cluster-steps] 
    [--hmm-identity] [--hmm-coverage] [--hmm-evalue] [--hmm-sensitivity] [--hmm-cluster-mode] [--hmm-cluster-steps] 
    [--skip-hmm] [-v] [-d] [-a] [-p] [-c CPUS] [-o OUTDIR] infile [infile ...]

    Assort phage protein sequences into phamilies of homologs using MMseqs2.

    positional arguments:
      infile                path to input file(s) in FASTA or Genbank flatfile format

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         print progress messages to the console
      -d, --debug           run in debug mode
      -a, --align-phams     use Clustal Omega to align phams (could take awhile...)
      -p, --pangenome       pangenome analysis à la Roary (only meaningful if given one input file per genome)
      -c , --cpus           number of threads to use [default: 4]
      -o , --outdir         path to directory where output files should go

    MMseqs2 sequence-sequence clustering arguments:
      --identity            percent identity for sequence-sequence clustering [default: 35.0%]
      --coverage            percent coverage for sequence-sequence clustering [default: 80.0%]
      --evalue              E-value threshold for sequence-sequence clustering [default: 0.001]
      --sensitivity         sensitivity: 1 favors speed, 7 favors sensitivity [default: 7]
      --cluster-mode        clustering algorithm [default: 0]
      --cluster-steps       number of steps for sequence-sequence clustering to proceed in [default: 1]

    MMseqs2 profile-sequence clustering arguments:
      --hmm-identity        percent identity for profile-consensus clustering [default: 15.0%]
      --hmm-coverage        percent coverage for profile-consensus clustering [default: 70.0%]
      --hmm-evalue          E-value threshold for profile-consensus clustering [default: 0.001]
      --hmm-sensitivity     sensitivity: 1 favors speed, 7 favors sensitivity [default: 7]
      --hmm-cluster-mode    clustering algorithm [default: 0]
      --hmm-cluster-steps   number of steps for profile-consensus clustering to proceed in [default: 3]
      --skip-hmm            do not perform profile-consensus clustering

    Steinegger M. and Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data 
    sets. Nature Biotechnology, 2017. doi: 10.1038/nbt.3988

The only required argument is the path to a single multiple-FASTA file, for example:

    phammseqs my_genes.faa

This will perform pham assembly, and write each resultant pham to a FASTA file found in a new folder in the working 
directory called `phammseqs__[day]_[month]_[year]` (this will resolve to the date on which pham assembly was performed).

An alternate output path can be specified with the `-o` argument:

    phammseqs my_genes.faa -o ~/Desktop/phammseqs_results

This will do the same as before, the output files will be now found in `~/Desktop/phammseqs_results` rather than the
directory the program was invoked from.

If your dataset is a pangenome or metagenome with many FASTA files (e.g. one file per genome), you can specify multiple
input files by simply putting their paths one after the next:

    phammseqs genome1.faa genome2.faa genome3.faa ... genomeN.faa -o ~/Desktop/phammseqs_results

or if all these genomes are in the same directory:

    phammseqs /path/to/genome/fastas/*.faa -o ~/Desktop/phammseqs_results
    
Each input file is treated separately, so you can even mix FASTA and Genbank flatfiles in the same run:

    phammseqs /path/to/genome/fastas/*.faa /path/to/genome/genbanks/*.gbk -o ~/Desktop/phammseqs_results

If you want to produce a multiple sequence alignment for each pham, the phammseqs program can accomplish this using 
a local copy of the program `clustalo` - simply use the `-a`/`--align-phams` argument:

    phammseqs my_genes.faa -o ~/Desktop/phammseqs_results -a -v

The `-v` argument will make the program print progress messages to the console as it runs, for example:

    Parsing protein sequences from input files...
    Found 404954 translations in 1 file(s)...
    Creating MMseqs2 database...
    Performing sequence-sequence clustering...
    Parsing sequence-sequence phams...
    Building profiles from sequence-sequence phams...
    Extracting consensus sequences from profiles...
    Performing profile-consensus clustering...
    Parsing profile-consensus phams...
    Found 27358 phamilies in dataset...
    Computing phamily alignments with Clustal Omega...
    [############                                     ] 25%

This may be especially helpful on very large or highly diverse datasets.

# Using PhaMMseqs as a library

For most simple use cases, the following import statement should suffice.

    from phammseqs import *

This will import two classes and two high-level functions into the namespace: `SequenceDB`, `Pham`, `assemble_phams`, 
and `merge_seq_hmm_phams`. Initialize a `SequenceDB` to begin:

    db = SequenceDB()

Load the contents of a FASTA file into a SequenceDB instance like this:

    db.load("/path/to/file.fasta")

The FASTA parser does not care if your file is in 2-line or multi-line/wrapped FASTA format (or mixed). So long as the
file does not contain duplicate headers, the `load()` function should work. If your FASTA file contains duplicate 
headers, a `ValueError` will be raised, so this alternate database load strategy can be used:

    from phammseqs.fileio import read_fasta
    
    for header, sequence in read_fasta("/path/to/file.fasta"):
        try:
            db.add_gene(header, sequence)
        except ValueError as err:
            print(err)      # will print the error message but keep adding genes to db

The size of the database can be queried two different ways:

    len(db)                 # number of genes in the database
    len(db.translations)    # number of non-redundant genes in the database

The complete database can be iterated over like so:

    for geneid, translation in db:
        print(f">{geneid}\n{translation}\n")

Alternatively, just the non-redundant sequences can be iterated:

    for geneid, translation in db.nr_genes:
        print(f">{geneid}\n{translation}\n")

A FASTA file can be written containing just the non-redundant sequences like this:

    db.write("/path/to/file.fasta", nr=True)    # nr=False for all genes

Pham assembly is pretty easy with the `assemble_phams` function. You'll need to define the clustering parameters in two
dictionaries (one for sequence-sequence clustering, one for profile-consensus clustering):

    seq_params = {"identity": 35, "coverage": 80, "evalue": 0.001,
                  "sensitivity": 7, "cluster_mode": 0, "cluster_steps": 1}
    hmm_params = {"identity": 15, "coverage": 70, "evalue": 0.001,
                  "sensitivity": 7, "cluster_mode": 0, "cluster_steps": 3}
    
    phams = assemble_phams(db, seq_params, hmm_params)

If profile-consensus clustering is not needed, simply omit the hmm_params:

    phams = assemble_phams(db, seq_params)

`phams` is a list of `Pham` objects, which have a few useful attributes and functions in addition to all those found 
in the parent class (`SequenceDB`), for example:

    pham = phams[0]
    pham.minimum_length     # return length of shortest gene -> int
    pham.average_length     # return average gene length     -> float
    pham.maximum_length     # return length of longest gene  -> int
    pham.is_orpham          # is this pham an orpham?        -> bool

If one has reason to do so, phams can easily be merged:

    first, second = phams[0], phams[1]
    len(first), len(second) # 775, 682
    first.merge(second)     # copies every sequence from second into first

Far more advanced workflows can be built leveraging the other submodules in phammseqs. For example, one can directly 
invoke MMseqs2 pipelines, or run ClustalO to generate MSAs of individual FASTA files of interest.