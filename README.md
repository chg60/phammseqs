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
required to run it. However, if your system does not have Python 3.6 or higher, you will need to obtain it 
[here](https://www.python.org/downloads/).

Once all that is done, you can obtain the phamerate package from PyPI using pip:

    pip3 install phamerate

# Basic Usage
With all dependencies met, you can run phamerate by invoking it with the `-h` option (to print the help menu):

    phamerate -h

Which should print something like:

    usage: phamerate [-h] [--cluster-mode] [--sensitivity] [--identity] [--coverage] [--evalue] [--no-hmm] [--hmm-identity] [--hmm-coverage] [--hmm-evalue] [-c] [-v] [-o] [-t] [-a] infile [infile ...]
    
    Assort phage protein sequences into phamilies using MMseqs2.
    
    positional arguments:
      infile             path to input file(s) in FASTA format
    
    optional arguments:
      -h, --help         show this help message and exit
      -c , --cpus        number of threads to use [default: 8]
      -v, --verbose      print progress messages
      -o , --outdir      path to directory where output files should go [default: /Users/your_username]
      -t , --tmpdir      path where temporary file I/O should occur [default: /tmp/phamerate]
      -a, --align-phams  use Clustal Omega to align phams (this could take awhile...)
    
    mmseqs arguments:
      --cluster-mode     clustering algorithm [default: 0]
      --sensitivity      sensitivity: 1.0 favors speed, 7.5 favors sensitivity [default: 4.0]
      --identity         percent identity for sequence-sequence clustering [default: 0.3]
      --coverage         percent coverage for sequence-sequence clustering [default: 0.85]
      --evalue           E-value threshold for sequence-sequence clustering [default: 0.001]
      --no-hmm           skip HMM clustering
      --hmm-identity     percent identity for consensus-HMM clustering [default: 0.25]
      --hmm-coverage     percent coverage for consensus-HMM clustering [default: 0.5]
      --hmm-evalue       E-value threshold for consensus-HMM clustering [default: 0.001]
    
    Steinegger M. and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, 2017. doi: 10.1038/nbt.3988

The only required argument is the path to a single multiple-FASTA file, for example:

    phamerate my_genes.faa

This will perform pham assembly, and create a directory `phamily_fastas` containing a multiple-FASTA file for each gene 
phamily in the input gene set.

An alternate output path can be specified with the `-o` argument:

    phamerate my_genes.faa -o ~/Desktop/phamerate_results

This will do the same as before, except `phamily_fastas` will be found in `~/Desktop/phamerate_results` rather than the
directory the program was invoked from.

If your dataset is a pangenome or metagenome with many FASTA files (e.g. one file per genome), you can specify multiple
input files by simply putting their paths one after the next:

    phamerate genome1.faa genome2.faa genome3.faa ... genomeN.faa -o ~/Desktop/phamerate_results

or if all these genomes are in the same directory:

    phamerate /path/to/genome/fastas/*.faa -o ~/Desktop/phamerate_results

If you want to have an MSA for each pham, the phamerate program can accomplish this using `clustalo` - simply use the 
`-a` argument:

    phamerate my_genes.faa -o ~/Desktop/phamerate_results -a -v

The `-v` argument will make the program print progress messages to the console as it runs, for example:

    Found 378159 translations in 1 file(s)...
    Creating MMseqs2 database...
    Performing sequence-sequence clustering...
    Parsing first iteration phams...
    Building HMMs from pre-phams...
    Extracting consensus sequences from HMMs...
    Performing consensus-HMM clustering...
    Parsing second iteration phams...
    Found 22897 phamilies in dataset...
    Aligning phams with Clustal Omega...
    [############                                     ] 25%

This may be especially helpful on large datasets, as the progressbar updates to show what fraction of alignments have 
been computed. This should give you a sense of whether you have time to make a cup of coffee while it finishes...

# Advanced Usage
For folks with large bacterial pan-genomes to analyze, you may find that the BLAST-based method used by 
[Roary](https://github.com/sanger-pathogens/Roary) is too slow for your needs. In this case, phamerate may be able to 
help, by raising the `--identity` threshold to 0.95 (same as the 95% identity threshold used by Roary) and supplying 
the `--no-hmm` argument, as you won't be searching for remote homologs:

    phamerate my_genes.faa -o ~/Desktop/phamerate_results -a -v --identity 0.95 --no-hmm

Which will print something like:

    Found 378159 translations in 3829 files...
    Creating MMseqs2 database...
    Performing sequence-sequence clustering...
    Parsing first iteration phams...
    Found 93674 phamilies in dataset...
    Computing phamily alignments with Clustal Omega...
    [############                                     ] 25%

Fair warning: at present, phamerate does NOT make any effort to split paralogs out of gene phamilies.

# Future Releases
We would like to do the following in future releases:
* remove paralogs from phamilies                                        [enhancement]
* export figure showing marginal pan-genome (each new genome adds...)   [enhancement]
* create tree(s) based on well-conserved genes                          [enhancement]