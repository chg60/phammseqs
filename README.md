# phamerate

The phamerate package facilitates pham assembly using [MMseqs2](https://www.mmseqs.com). Default parameters have
been carefully tuned for rapid, accurate exploration of the bacteriophage protein sequence space.

# Conda installation

The easiest way to install the phamerate package and its dependencies is through the Anaconda/Miniconda package manager: 

    conda create -n phamerate-env python=3.9 -y && conda activate phamerate-env
    conda install -c bioconda -c conda-forge mmseqs2=13.45111 clustalo -y
    pip3 install phamerate 

# Manual installation

If you don't have some flavor of conda available (and don't want to install it...) you may follow the instructions
[here](https://github.com/soedinglab/mmseqs2#installation) to manually install `mmseqs`. An optional dependency,
`clustalo` can be manually installed following the instructions [here](http://www.clustal.org/omega/). 
Most modern operating systems also ship with Python3, the programming language used to develop this package, and 
required to run_clustalo it. However, if your system does not have Python 3.7 or higher, you will need to obtain it 
[here](https://www.python.org/downloads/).

Once all that is done, you can obtain the phamerate package from PyPI using pip:

    pip3 install phamerate

# Basic Usage

With all dependencies met, you can run phamerate by invoking it with the `-h` option (to print the help menu):

    phamerate -h

Which should print something like:

    usage: phamerate [-h] [--identity] [--coverage] [--evalue] [--sensitivity] [--cluster-mode] [--cluster-steps] 
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

    phamerate my_genes.faa

This will perform pham assembly, and write each resultant pham to a FASTA file found in a new folder in the working 
directory called `phamerate__[day]_[month]_[year]` (this will resolve to the date on which pham assembly was performed).

An alternate output path can be specified with the `-o` argument:

    phamerate my_genes.faa -o ~/Desktop/phamerate_results

This will do the same as before, the output files will be now found in `~/Desktop/phamerate_results` rather than the
directory the program was invoked from.

If your dataset is a pangenome or metagenome with many FASTA files (e.g. one file per genome), you can specify multiple
input files by simply putting their paths one after the next:

    phamerate genome1.faa genome2.faa genome3.faa ... genomeN.faa -o ~/Desktop/phamerate_results

or if all these genomes are in the same directory:

    phamerate /path/to/genome/fastas/*.faa -o ~/Desktop/phamerate_results
    
Each input file is treated separately, so you can even mix FASTA and Genbank flatfiles in the same run!

    phamerate /path/to/genome/fastas/*.faa /path/to/genome/genbanks/*.gbk -o ~/Desktop/phamerate_results

If you want to produce a multiple sequence alignment for each pham, the phamerate program can accomplish this using 
a local copy of the program `clustalo` - simply use the `-a`/`--align-phams` argument:

    phamerate my_genes.faa -o ~/Desktop/phamerate_results -a -v

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

# Can I use it as a stand-in for my favorite pan-genome analysis tool?

For folks with large bacterial pan-genomes to analyze, you may find that the BLAST-based method used by 
[Roary](https://github.com/sanger-pathogens/Roary) or other such tools is too slow for your needs; in this case, 
phamerate may be able to help. By raising the `--identity` threshold (for example to the 95% identity threshold used 
by Roary) and supplying the `--skip-hmm` argument (no searching for remote homologs) and the `-p` argument to generate 
some tabular output files very similar to those produced by Roary:

    phamerate ~/Desktop/Actino_phages/*.faa -o ~/Desktop/phamerate_results -a -v --identity 95 --skip-hmm -p

Which as it runs and completes steps should print something like:

    Parsing protein sequences from input files...
    Found XXXXXX translations in YYY files...
    Creating MMseqs2 database...
    Performing sequence-sequence clustering...
    Parsing sequence-sequence phams...
    Found ZZZZ phamilies in dataset...
    Computing phamily alignments with Clustal Omega...
    [##################################################] 100%
    Performing basic pan/meta-genome analyses...
    Done!

Be forewarned that at present, phamerate does NOT make any effort to split paralogs out of gene phamilies, so if that 
is something that matters for your analyses you'll need to find another way to split over-clustered phams.
