"""Functions for performing pham assembly with MMseqs2."""

from phamerate import subprocess

VERSION = "13.45111"


class MMseqs2Error(Exception):
    """Error caused by running MMseqs2."""
    pass


def get_mmseqs_version():
    """Return the MMseqs2 version string.

    :return: version
    """
    command = "mmseqs version"

    try:
        stdout, stderr = subprocess.run(command)
        if stderr:
            print(stderr)
            raise MMseqs2Error(f"command failed: {command}")
        return stdout.rstrip()
    except FileNotFoundError:
        return None


def createdb(fasta, seq_db, debug=False):
    """Run 'mmseqs createdb' command.

    :param fasta: path to the FASTA file to convert
    :type fasta: pathlib.Path
    :param seq_db: path to the sequence database to create
    :type seq_db: pathlib.Path
    :param debug: print verbose output to console
    :type debug: bool
    """
    command = f"mmseqs createdb {fasta} {seq_db} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return seq_db


def cluster(seq_db, clu_db, tmp_dir, identity, coverage, e_value, sens,
            max_seqs=1000, cluster_mode=0, cluster_reassign=True,
            cluster_steps=1, alignment_mode=3, coverage_mode=0, cpus=1,
            debug=False):
    """Run 'mmseqs cluster' command.

    Arguments are mapped to the following MMseqs2 arguments:

    * identity = `--min-seq-id`
    * coverage = `-c`
    * e_value = `-e`
    * sens = `-s`
    * max_seqs = `--max-seqs`
    * cluster_mode = `--cluster-mode`
    * cluster_reassign = `--cluster-reassign`
    * cluster_steps = `--cluster-steps`
    * alignment_mode = `--alignment-mode`
    * coverage_mode = `--cov-mode`
    * cpus = `--threads`

    :param seq_db: path to the sequence database to cluster
    :type seq_db: pathlib.Path
    :param clu_db: path to write_fasta the clustered sequence database
    :type clu_db: pathlib.Path
    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param identity: percent identity threshold for clustering
    :type identity: float
    :param coverage: percent coverage threshold for clustering
    :type coverage: float
    :param e_value: significance threshold for clustering
    :type e_value: float
    :param sens: sensitivity parameter for seeding alignments
    :type sens: float
    :param max_seqs: number of hits to store per query
    :type max_seqs: int
    :param cluster_mode: clustering algorithm to use
    :type cluster_mode: int
    :param cluster_reassign: reassign cluster members that violate
        thresholds
    :type cluster_reassign: bool
    :param cluster_steps: number of clustering steps to proceed with
    :type cluster_steps: int
    :param alignment_mode: what kind of alignments to generate
    :type alignment_mode: int
    :param coverage_mode: coverage calculation to use
    :type coverage_mode: int
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to console
    :type debug: bool
    """
    if identity > 1:
        identity /= 100.0
    if coverage > 1:
        coverage /= 100.0

    command = f"mmseqs cluster {seq_db} {clu_db} {tmp_dir} " \
              f"--min-seq-id {identity} -c {coverage} -e {e_value} " \
              f"-s {sens} --max-seqs {max_seqs} --cluster-mode " \
              f"{cluster_mode} --cluster-steps {cluster_steps} " \
              f"--alignment-mode {alignment_mode} --cov-mode " \
              f"{coverage_mode}  --threads {cpus} -v 3 "

    if cluster_reassign:
        command += " --cluster-reassign"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return clu_db


def result2profile(seq_db, clu_db, pro_db, cpus=1, debug=False):
    """Run 'mmseqs result2profile' command.

    :param seq_db: path to the sequence database to cluster
    :type seq_db: pathlib.Path
    :param clu_db: path to the clustered sequence database
    :type clu_db: pathlib.Path
    :param pro_db: path to the desired profile database
    :type pro_db: pathlib.Path
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to console
    :type debug: bool
    """
    command = f"mmseqs result2profile {seq_db} {seq_db} {clu_db} {pro_db} " \
              f"--threads {cpus} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return pro_db


def profile2consensus(pro_db, con_db, cpus=1, debug=False):
    """Run 'mmseqs profile2consensus' command.

    :param pro_db: path to a profile database
    :type pro_db: pathlib.Path
    :param con_db: path to the desired consensus sequence database
    :type con_db: pathlib.Path
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to console
    :type debug: bool
    """
    command = f"mmseqs profile2consensus {pro_db} {con_db} " \
              f"--threads {cpus} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return con_db


def search(pro_db, con_db, aln_db, tmp_dir, identity, coverage, e_value, sens,
           max_seqs=1000, cluster_steps=1, alignment_mode=3, coverage_mode=0,
           cpus=1, debug=False):
    """Run 'mmseqs search' command.

    Arguments are mapped to the following MMseqs2 arguments:

    * identity = `--min-seq-id`
    * coverage = `-c`, `--cov`
    * e_value = `-e`, `--e-profile`
    * sens = `-s`
    * max_seqs = `--max-seqs`
    * cluster_steps = `--cluster-steps`
    * alignment_mode = `--alignment-mode`
    * coverage_mode = `--cov-mode`
    * cpus = `--threads`

    :param pro_db: path to a profile database
    :type pro_db: pathlib.Path
    :param con_db: path to a consensus sequence database
    :type con_db: pathlib.Path
    :param aln_db: path to the desired alignment database
    :type aln_db: pathlib.Path
    :param tmp_dir: directory MMseqs2 can use as scratch space
    :type tmp_dir: pathlib.Path
    :param identity: percent identity threshold
    :type identity: float
    :param coverage: percent coverage threshold
    :type coverage: float
    :param e_value: significance threshold
    :type e_value: float
    :param sens: sensitivity parameter for seeding alignments
    :type sens: float
    :param max_seqs: maximum number of hits to store per query
    :type max_seqs: int
    :param cluster_steps: num of iterative profile-consensus searches
        to perform
    :type cluster_steps: int
    :param alignment_mode: MMseqs2 alignment mode
    :type alignment_mode: int
    :param coverage_mode: MMseqs2 coverage mode
    :type coverage_mode: int
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to console
    :type debug: bool
    """
    if identity > 1:
        identity /= 100.0
    if coverage > 1:
        coverage /= 100.0

    command = f"mmseqs search {pro_db} {con_db} {aln_db} {tmp_dir} " \
              f"--min-seq-id {identity} -c {coverage} --cov {coverage} " \
              f"-e {e_value} --e-profile {e_value} -s {sens} --max-seqs " \
              f"{max_seqs} --num-iterations {cluster_steps} " \
              f"--alignment-mode {alignment_mode} --cov-mode " \
              f"{coverage_mode} --add-self-matches --threads {cpus} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return aln_db


def clust(con_db, aln_db, res_db, cluster_mode, cpus=1, debug=False):
    """Run 'mmseqs clust' command.

    :param con_db: directory MMseqs2 can use as scratch space
    :type con_db: pathlib.Path
    :param aln_db: directory MMseqs2 can use as scratch space
    :type aln_db: pathlib.Path
    :param res_db: directory MMseqs2 can use as scratch space
    :type res_db: pathlib.Path
    :param cluster_mode: clustering algorithm to use
    :type cluster_mode: int
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to console
    :type debug: bool
    """
    command = f"mmseqs clust {con_db} {aln_db} {res_db} " \
              f"--cluster-mode {cluster_mode} --threads {cpus} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return res_db


def createseqfiledb(seq_db, clu_db, sf_db, cpus=1, debug=False):
    """Run 'mmseqs createseqfiledb' command.

    :param seq_db: path to the sequence database
    :type seq_db: pathlib.Path
    :param clu_db: path to the clustered sequence database
    :type clu_db: pathlib.Path
    :param sf_db: path where the seqfiledb should be written
    :type sf_db: pathlib.Path
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param debug: print verbose output to the console
    :type debug: bool
    :return: sf_db
    """
    command = f"mmseqs createseqfiledb {seq_db} {clu_db} {sf_db} " \
              f"--threads {cpus} -v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return sf_db


def result2flat(query_db, target_db, sf_db, outfile, debug=False):
    """Run 'mmseqs result2flat' command after clustering.

    :param query_db: the query database used for clustering
    :type query_db: pathlib.Path
    :param target_db: the target database used for clustering
    :type target_db: pathlib.Path
    :param sf_db: the seqfiledb created after clustering
    :type sf_db: pathlib.Path
    :param outfile: path where the FASTA-like cluster output should go
    :type outfile: pathlib.Path
    :param debug: print verbose output to the console
    :type debug: bool
    :return: outfile
    """
    command = f"mmseqs result2flat {query_db} {target_db} {sf_db} {outfile} " \
              f"-v 3"

    stdout, stderr = subprocess.run(command)
    if debug and stdout:
        print(stdout)
    if stderr:
        print(stderr)
        raise MMseqs2Error(f"command failed: {command}")

    return outfile


def parse_output(outfile):
    """Parse the indicated MMseqs2 FASTA-like file into a dictionary of
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

        # Need to write_fasta the last pham into the dictionary
        phams[pham_name] = pham_translations

    # Pham -1 is a placeholder
    phams.pop(-1)

    return phams
