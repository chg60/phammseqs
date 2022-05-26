"""Assort phage protein sequences into phamilies of homologs using MMseqs2."""

import pathlib
import shutil
import sys
import tempfile

from phamerate.SequenceDB import Pham, SequenceDB
from phamerate import cli, clustalo, fileio, mmseqs
from phamerate import multiprocess as mp
from phamerate.pangenome import analyze_pangenome

# Check MMseqs2 version - exit if not found or wrong version
mmseqs_version = mmseqs.get_mmseqs_version()
if not mmseqs_version:
    print(f"\nMMseqs2 not found - please install version "
          f"{mmseqs.VERSION} and try again\n")
    sys.exit(1)
elif mmseqs_version != mmseqs.VERSION:
    print(f"\nFound wrong version of MMseqs2")
    print(f"Installed version: {mmseqs_version}")
    print(f"Required version:  {mmseqs.VERSION}\n")
    sys.exit(1)


def create_database(input_files):
    """Parse the input file(s) into a SequenceDB object.

    :param input_files: the input file(s) to parse sequences from
    :type input_files: list[pathlib.Path]
    :return: database
    """
    database = SequenceDB()
    for input_file in input_files:
        fmt = fileio.sniff_format(input_file)
        if fmt == "fasta":
            sequences = fileio.read_fasta(input_file)
        elif fmt == "genbank":
            sequences = fileio.read_genbank(input_file)
        else:
            print(f"Unknown file format for {str(input_file)} - skipping")
            continue

        for geneid, translation in sequences:
            try:
                database.add_gene(geneid, translation.upper())
            except ValueError as err:
                print(err)

    return database


def merge_seq_hmm_phams(hmm_phams, seq_phams, database):
    """Merge sequence-sequence and profile-consensus phams.

    :param hmm_phams: profile-consensus clusters
    :type hmm_phams: dict
    :param seq_phams: sequence-sequence clusters
    :type seq_phams: dict
    :param database: the database whose sequences were clustered
    :type database: SequenceDB
    :return: merged_phams
    """
    phams, lookup = list(), dict()
    for key, translations in seq_phams.items():
        pham = Pham()
        for translation in translations:
            lookup[translation] = key
            for geneid in database.get_translation_geneids(translation):
                pham.add_gene(geneid, translation)
        phams.append(pham)

    if not hmm_phams:
        return sorted(phams, reverse=True)

    merged_phams = list()
    for key, translations in hmm_phams.items():
        pham = Pham()
        for translation in translations:
            pre_pham_key = lookup[translation]
            pham.merge(phams[pre_pham_key])
        merged_phams.append(pham)

    return sorted(merged_phams, reverse=True)


def assemble_phams(db, seq_params, hmm_params=None, cpus=1,
                   verbose=False, debug=False):
    """Perform pham assembly on the indicated database.

    If `hmm_params` is `None`, only sequence-sequence clustering will be
    done.

    :param db: the database whose sequences need to be assorted into
        gene phamilies
    :type db: SequenceDB
    :param seq_params: sequence-sequence clustering parameters
    :type seq_params: dict
    :param hmm_params: profile-consensus clustering parameters
    :type hmm_params: dict
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param verbose: print verbose output to console
    :type verbose: bool
    :param debug: run_clustalo MMseqs2 commands in debug mode
    :type debug: bool
    :return: phams
    """
    # Set up temporary directory and name our temporary files
    tmp_dir = pathlib.Path(tempfile.mkdtemp(prefix="phamerate-"))

    nr_fasta = tmp_dir.joinpath("nr_genes.fasta")

    seq_db = tmp_dir.joinpath("sequenceDB")
    clu_db = tmp_dir.joinpath("seqClusterDB")
    sf_db1 = tmp_dir.joinpath("seqSeqfileDB")
    seq_out = tmp_dir.joinpath("seqPhams.txt")

    pro_db = tmp_dir.joinpath("profileDB")
    con_db = tmp_dir.joinpath("consensusDB")
    aln_db = tmp_dir.joinpath("alignDB")
    res_db = tmp_dir.joinpath("hmmClusterDB")
    sf_db2 = tmp_dir.joinpath("hmmSeqfileDB")
    hmm_out = tmp_dir.joinpath("hmmPhams.txt")

    try:
        if verbose:
            print("Creating MMseqs2 database...")

        db.write(nr_fasta, nr_only=True)
        mmseqs.createdb(nr_fasta, seq_db, debug=debug)

        if verbose:
            print("Performing sequence-sequence clustering...")

        # Unpack sequence-sequence clustering parameters, then run_clustalo
        i, c = seq_params["identity"], seq_params["coverage"]
        e, s = seq_params["evalue"], seq_params["sensitivity"]
        m, n = seq_params["cluster_mode"], seq_params["cluster_steps"]
        mmseqs.cluster(seq_db, clu_db, tmp_dir, identity=i, coverage=c,
                       e_value=e, sens=s, cluster_mode=m, cluster_steps=n,
                       cpus=cpus, debug=debug)

        if verbose:
            print("Parsing sequence-sequence phams...")

        mmseqs.createseqfiledb(seq_db, clu_db, sf_db1, cpus=cpus,
                               debug=debug)
        mmseqs.result2flat(seq_db, seq_db, sf_db1, seq_out, debug=debug)

        seq_phams = mmseqs.parse_output(seq_out)

        if hmm_params:
            if verbose:
                print("Building profiles from sequence-sequence phams...")

            mmseqs.result2profile(seq_db, clu_db, pro_db, cpus=cpus,
                                  debug=debug)
            if verbose:
                print("Extracting consensus sequences from profiles...")

            mmseqs.profile2consensus(pro_db, con_db, cpus=cpus,
                                     debug=debug)

            if verbose:
                print("Performing profile-consensus clustering...")

            i, c = hmm_params["identity"], hmm_params["coverage"]
            e, s = hmm_params["evalue"], hmm_params["sensitivity"]
            m, n = hmm_params["cluster_mode"], hmm_params["cluster_steps"]
            mmseqs.search(pro_db, con_db, aln_db, tmp_dir, identity=i,
                          coverage=c, e_value=e, sens=s, cluster_steps=n,
                          cpus=cpus, debug=debug)

            mmseqs.clust(con_db, aln_db, res_db, m, cpus=cpus, debug=debug)

            if verbose:
                print("Parsing profile-consensus phams...")

            mmseqs.createseqfiledb(seq_db, res_db, sf_db2, cpus=cpus,
                                   debug=debug)
            mmseqs.result2flat(pro_db, con_db, sf_db2, hmm_out,
                               debug=debug)

            hmm_phams = mmseqs.parse_output(hmm_out)
        else:
            hmm_phams = None
    except mmseqs.MMseqs2Error as err:
        # Catch MMseqs2 errors, so we can copy temporary files to ~/
        if debug:
            dst = pathlib.Path().home().joinpath(f"{tmp_dir.name}")
            shutil.copytree(src=tmp_dir, dst=dst)
            print(f"Temporary files are at {dst}")
        # Re-raise the error
        raise err
    finally:
        # Always clear temporary directory regardless of error status
        shutil.rmtree(tmp_dir)

    return merge_seq_hmm_phams(hmm_phams, seq_phams, db)


def write_phams(phams, fasta_dir, nr_only=False):
    """Dump each pham to a FASTA file. Return the path and effective
    size of each pham.

    :param phams: the assembled phams
    :type phams: list[Pham]
    :param fasta_dir: the directory where files should be written
    :type fasta_dir: pathlib.Path
    :param nr_only: only write_fasta non-redundant sequences to file
    :type nr_only: bool
    :return: pham_fastas
    """
    pham_fastas = list()
    for i, pham in enumerate(phams):
        fasta_path = fasta_dir.joinpath(f"pham_{i + 1}.faa")
        pham.write(fasta_path, nr_only=nr_only)
        pham_fastas.append(fasta_path)

    return pham_fastas


def align_phams(phams, fasta_dir, align_dir, cpus=1, verbose=False):
    """Generate pham MSAs with Clustal Omega.

    :param phams: the phams to align
    :type phams: list[Pham]
    :param fasta_dir: path where the fastas should go
    :type fasta_dir: pathlib.Path
    :param align_dir: path where the alignments should go
    :type align_dir: pathlib.Path
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param verbose: print progress bar
    :type verbose: bool
    """
    fasta_paths, align_paths = list(), list()
    jobs = list()
    for i, pham in enumerate(phams):
        fasta_path = fasta_dir.joinpath(f"pham_{i + 1}.faa")
        fasta_paths.append(fasta_path)

        # Not sensible to align phams with only 1 member
        if len(pham) == 1:
            pham.write(fasta_path)
            continue

        align_path = align_dir.joinpath(f"pham_{i + 1}.aln")
        align_paths.append(align_path)

        # If all sequences are the same, they are already aligned
        if pham.n_eff == 1:
            pham.write(fasta_path)
            pham.write(align_path)
            continue

        # All other cases we need to align
        pham.write(fasta_path, nr_only=True)
        jobs.append((fasta_path, align_path))

    results = mp.parallelize(jobs, cpus, clustalo.run_clustalo, verbose)
    for fasta_path, align_path in results:
        i = int(align_path.stem.split("_")[-1]) - 1
        pham = phams[i]

        pham.update_msa(align_path)
        pham.write(fasta_path)

    return fasta_paths, align_paths


def phamerate(infiles, seq_params, hmm_params, cpus=1, verbose=False,
              debug=False):
    """Perform pham assembly with the sequences found in `infiles` and
    parameters specified in `seq_params` and `hmm_params` and results
    written to `outdir`.

    :param infiles: input files containing sequences to cluster
    :type infiles: list[pathlib.Path]
    :param seq_params: sequence-sequence clustering parameters
    :type seq_params: dict
    :param hmm_params: profile-consensus clustering parameters
    :type hmm_params: dict
    :param cpus: number of CPU cores to use
    :type cpus: int
    :param verbose: print progress messages to the console
    :type verbose: bool
    :param debug: print debug messages to the console
    :type debug: bool
    """
    if verbose:
        print("Parsing protein sequences from input files...")

    database = create_database(infiles)

    if verbose:
        print(f"Found {len(database)} translations in {len(infiles)} files...")

    phams = assemble_phams(db=database, seq_params=seq_params,
                           hmm_params=hmm_params, cpus=cpus,
                           verbose=verbose, debug=debug)

    return phams


def main():
    """Commandline entry point for this module."""
    # Invoke the help menu if phamerate was run_clustalo without args
    if len(sys.argv) == 1:
        sys.argv.append("-h")
    args = cli.parse_args()

    infiles = args.infile
    outdir = args.outdir
    if not outdir.is_dir():
        outdir.mkdir(parents=True)

    skip_hmm = args.skip_hmm
    debug = args.debug
    verbose = args.verbose or debug
    cpus = args.cpus
    do_align = args.align_phams
    do_pangenome = args.pangenome

    # If user specified align, check for compatible Clustal Omega version
    if do_align:
        # Always warn regardless of verbose flag
        clustalo_version = clustalo.get_clustalo_version()
        if not clustalo_version:
            print(f"\nClustal Omega not found - pham alignments disabled...\n")
            do_align = False
        elif clustalo_version != clustalo.VERSION:
            print(f"\nFound wrong version of Clustal Omega")
            print(f"Installed version: {clustalo_version}")
            print(f"Required version: {clustalo.VERSION}")
            print(f"Pham alignments disabled...\n")
            do_align = False

    # Set up MMseqs2 parameters
    seq_params = {"identity":         args.identity,
                  "coverage":         args.coverage,
                  "evalue":           args.evalue,
                  "sensitivity":      args.sensitivity,
                  "cluster_mode":     args.cluster_mode,
                  "cluster_steps":    args.cluster_steps}

    if skip_hmm:
        hmm_params = None
    else:
        hmm_params = {"identity":      args.hmm_identity,
                      "coverage":      args.hmm_coverage,
                      "evalue":        args.hmm_evalue,
                      "sensitivity":   args.sensitivity,
                      "cluster_mode":  args.hmm_cluster_mode,
                      "cluster_steps": args.hmm_cluster_steps}

    phams = phamerate(infiles, seq_params, hmm_params, cpus=cpus,
                      verbose=verbose, debug=debug)

    if verbose:
        print(f"Found {len(phams)} phamilies in dataset...")

    fasta_dir = outdir.joinpath("pham_fastas")
    if not fasta_dir.is_dir():
        fasta_dir.mkdir()

    if not do_align:
        write_phams(phams, fasta_dir)
    else:
        if verbose:
            print(f"Computing phamily alignments with Clustal Omega...")

        align_dir = outdir.joinpath("pham_aligns")
        if not align_dir.is_dir():
            align_dir.mkdir()

        align_phams(phams, fasta_dir, align_dir, cpus=cpus, verbose=verbose)

    if do_pangenome:
        if len(infiles) == 1:
            print("Cannot perform pangenome analysis on single input file")
        elif len(infiles) > 1:
            if verbose:
                print("Analyzing pangenome...")

            analyze_pangenome(phams, infiles, outdir)

    print("Done!")


if __name__ == "__main__":
    main()
