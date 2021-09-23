"""High-level utility functions for pham assembly and post-processing."""

from phamerate.classes import SequenceDB, Pham
from phamerate.fasta import parse_fasta, write_fasta
from phamerate.mmseqs import *


def create_database(input_files):
    """
    Walk a list of FASTA files and map non-redundant translations
    to the genomes and labels associated with them.

    :param input_files: FASTA file(s) to map translations from
    :type input_files: list of pathlib.Path
    :return: dataset
    """
    database = SequenceDB()
    for input_file in input_files:
        genome = input_file.name
        headers, sequences = parse_fasta(input_file)
        for geneid, translation in zip(headers, sequences):
            database.add_gene(geneid, genome, translation)

    return database


def dump_nr_translations(database, filepath):
    """
    Dump the non-redundant translations to a FASTA file for clustering.

    :param database: a database containing genes
    :type database: SequenceDB
    :param filepath: file to write non-redundant translations to
    :type filepath: pathlib.Path
    """
    headers, sequences = list(), list()
    for translation, data in database:
        headers.append(data[0][1])
        sequences.append(translation)

    write_fasta(headers, sequences, filepath)


def merge_pre_and_hmm_phams(hmm_phams, pre_phams, database):
    """
    Merge first and second iteration phams.

    :param hmm_phams: clustered consensus sequences
    :type hmm_phams: dict
    :param pre_phams: clustered sequences (used to generate hmms)
    :type pre_phams: dict
    :param database:
    :type database: SequenceDB
    :return: merged_phams
    """
    phams, lookup = list(), dict()
    for key, translations in pre_phams.items():
        pham = Pham()
        for translation in translations:
            lookup[translation] = key
            data = database.get_data(translation)
            for genome, geneid in data:
                pham.add_gene(geneid, genome, translation)
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


def assemble_phams(database, first_iter, second_iter, threads,
                   verbose, skip_hmm, tmp_dir):
    """
    Perform pham assembly and return the resultant phams.

    :param database: the sequence database
    :type database: SequenceDB
    :param first_iter: sequence-sequence clustering parameters
    :param first_iter: tuple of float
    :param second_iter: consensus-HMM clustering parameters
    :type second_iter: tuple of float
    :param threads: number of threads to use
    :type threads: int
    :param verbose: print messages along the way?
    :type verbose: bool
    :param skip_hmm: skip HMM clustering?
    :type skip_hmm: bool
    :param tmp_dir: temporary directory
    :type tmp_dir: pathlib.Path
    :return: phams
    """
    if verbose:
        print("Creating MMseqs2 database...")

    nr_fasta = tmp_dir.joinpath("all_genes.fasta")
    iter_1_out = tmp_dir.joinpath("firstIter.txt")
    iter_2_out = tmp_dir.joinpath("secondIter.txt")

    dump_nr_translations(database, nr_fasta)
    mmseqs_createdb(nr_fasta, tmp_dir)

    if verbose:
        print("Performing sequence-sequence clustering...")

    clu_mode, sens, i, c, e = first_iter
    mmseqs_cluster(tmp_dir, clu_mode, sens, i, c, e, threads=threads)

    if verbose:
        print("Parsing first iteration phams...")

    mmseqs_first_iter_cleanup(tmp_dir, iter_1_out, threads)

    first_iter_phams = parse_mmseqs_output(iter_1_out)

    if not skip_hmm:
        if verbose:
            print("Building HMMs from pre-phams...")

        mmseqs_result2profile(tmp_dir, threads)

        if verbose:
            print("Extracting consensus sequences from HMMs...")

        mmseqs_profile2consensus(tmp_dir, threads)

        if verbose:
            print("Performing consensus-HMM clustering...")

        clu_mode, sens, i, c, e = second_iter
        mmseqs_search(tmp_dir, i, c, e, threads=threads)

        mmseqs_clust(tmp_dir, clu_mode, threads)

        if verbose:
            print("Parsing second iteration phams...")

        mmseqs_second_iter_cleanup(tmp_dir, iter_2_out, threads)

        second_iter_phams = parse_mmseqs_output(iter_2_out)
    else:
        second_iter_phams = None

    return merge_pre_and_hmm_phams(second_iter_phams, first_iter_phams,
                                   database)
