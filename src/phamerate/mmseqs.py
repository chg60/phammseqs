"""Functions for performing pham assembly with MMseqs2."""

from phamerate.cmd import run_command


def mmseqs_createdb(fasta, tmp_dir, debug=False):
    """
    Run 'mmseqs createdb' command.

    :param fasta: path to the FASTA file to convert
    :type fasta: pathlib.Path
    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs createdb {fasta} {tmp_dir}/sequenceDB -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)


def mmseqs_cluster(tmp_dir, cluster_mode, sens, identity, coverage, evalue,
                   max_seqs=1000, steps=1, align_mode=3, cov_mode=0, threads=1,
                   debug=False):
    """
    Run 'mmseqs cluster' command.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param identity: percent identity threshold for clustering
    :type identity: float
    :param coverage: percent coverage threshold for clustering
    :type coverage: float
    :param evalue: significance threshold for clustering
    :type evalue: float
    :param sens: sensitivity parameter for seeding alignments
    :type sens: float
    :param cluster_mode: clustering algorithm to use
    :type cluster_mode: int
    :param max_seqs: number of hits to store per query
    :type max_seqs: int
    :param steps: number of clustering steps to proceed with
    :type steps: int
    :param align_mode: what kind of alignments to generate
    :type align_mode: int
    :param cov_mode: coverage calculation to use
    :type cov_mode: int
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs cluster {tmp_dir}/sequenceDB {tmp_dir}/clusterDB " \
              f"{tmp_dir} --min-seq-id {identity} -c {coverage} -e {evalue} " \
              f"-s {sens} --max-seqs {max_seqs} --cluster-steps {steps} " \
              f"--threads {threads} --alignment-mode {align_mode} --cov-mode " \
              f"{cov_mode} --cluster-mode {cluster_mode} --cluster-reassign " \
              f"-v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_result2profile(tmp_dir, threads=1, debug=False):
    """
    Run 'mmseqs result2profile' command.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs result2profile {tmp_dir}/sequenceDB " \
              f"{tmp_dir}/sequenceDB {tmp_dir}/clusterDB {tmp_dir}/profileDB " \
              f"--threads {threads} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_profile2consensus(tmp_dir, threads=1, debug=False):
    """
    Run 'mmseqs profile2consensus' command.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs profile2consensus {tmp_dir}/profileDB " \
              f"{tmp_dir}/consensusDB --threads {threads} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_search(tmp_dir, identity, coverage, evalue, max_seqs=1000,
                  align_mode=3, cov_mode=0, threads=1, debug=False):
    """
    Return 'mmseqs search' command.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param threads: how many threads to use?
    :type threads: int
    :param identity: percent identity threshold
    :type identity: float
    :param coverage: percent coverage threshold
    :type coverage: float
    :param evalue: significance threshold
    :type evalue: float
    :param threads: number of threads to use
    :type threads: int
    :param max_seqs: maximum number of hits to store per query
    :type max_seqs: int
    :param align_mode: MMseqs2 alignment mode
    :type align_mode: int
    :param cov_mode: MMseqs2 coverage mode
    :type cov_mode: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs search {tmp_dir}/profileDB {tmp_dir}/consensusDB " \
              f"{tmp_dir}/alignDB {tmp_dir} --min-seq-id {identity} -c " \
              f"{coverage} --cov  {coverage} -e {evalue} --e-profile " \
              f"{evalue} --threads  {threads} --max-seqs {max_seqs} " \
              f"--alignment-mode {align_mode} --cov-mode {cov_mode} " \
              f"--add-self-matches -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_clust(tmp_dir, cluster_mode, threads=1, debug=False):
    """
    Return 'mmseqs clust' command.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param cluster_mode: which clustering algorithm to use?
    :type cluster_mode: int
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs clust {tmp_dir}/consensusDB {tmp_dir}/alignDB " \
              f"{tmp_dir}/resultDB --cluster-mode {cluster_mode} " \
              f"--threads {threads} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_first_iter_cleanup(tmp_dir, outfile, threads=1, debug=False):
    """
    Run 'mmseqs createseqfiledb' command after sequence-sequence
    clustering.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param outfile: file where phams should be written
    :type outfile: pathlib.Path
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs createseqfiledb {tmp_dir}/sequenceDB " \
              f"{tmp_dir}/clusterDB {tmp_dir}/firstIterDB --threads " \
              f"{threads} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

    command = f"mmseqs result2flat {tmp_dir}/sequenceDB {tmp_dir}/sequenceDB " \
              f"{tmp_dir}/firstIterDB {outfile} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

def mmseqs_second_iter_cleanup(tmp_dir, outfile, threads=1, debug=False):
    """
    Run 'mmseqs createseqfiledb' command after consensus-HMM
    clustering.

    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param outfile: file where phams should be written
    :type outfile: pathlib.Path
    :param threads: how many threads to use?
    :type threads: int
    :param debug: run in debug mode (verbose output to console)
    :type debug: bool
    """
    command = f"mmseqs createseqfiledb {tmp_dir}/sequenceDB " \
              f"{tmp_dir}/resultDB {tmp_dir}/secondIterDB --threads " \
              f"{threads} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)

    command = f"mmseqs result2flat {tmp_dir}/profileDB {tmp_dir}/consensusDB " \
              f"{tmp_dir}/secondIterDB {outfile} -v 3"
    out, err = run_command(command, debug)
    if debug:
        print(out), print(err)


def parse_mmseqs_output(outfile):
    """
    Parse the indicated MMseqs2 FASTA-like file into a dictionary of
    integer-named phams.

    :param outfile: FASTA-like parseable output
    :type outfile: pathlib.Path
    :return: phams
    :rtype: dict
    """
    phams = {}
    pham_translations = list()
    pham_name = -1

    with open(f"{outfile}", "r") as fh:
        prior = fh.readline()
        current = fh.readline()

        # While loop to iterate until EOF
        while current:
            # If current and previous lines are header lines, new pham block
            if current.startswith(">") and prior.startswith(">"):
                phams[pham_name] = pham_translations
                pham_name += 1
                pham_translations = []

            # If the current line is a translation line, keep growing pham
            if not current.startswith(">"):
                pham_translations.append(current.rstrip())

            prior, current = current, fh.readline()

        # Need to dump the last pham into the dictionary
        phams[pham_name] = pham_translations

    # Pham -1 is a placeholder
    phams.pop(-1)

    return phams
