"""High-level utility functions for pham assembly and post-processing."""

from phamerate.classes import SequenceDB, Pham
from phamerate.fileio import sniff, parse_genbank, parse_fasta, write_fasta
from phamerate.mmseqs import *

HEADER = ["Gene", "Non unique Gene name", "No. isolates", "No. sequences",
          "Avg sequences per isolate", "Genome Fragment",
          "Order within Fragment", "Accessory Fragment",
          "Accessory Order with Fragment", "QC", "Min group size nuc",
          "Max group size nuc", "Avg group size nuc"]


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
        genome = input_file.stem
        fmt = sniff(input_file)
        if fmt == "fasta":
            headers, sequences = parse_fasta(input_file)
        elif fmt == "genbank":
            headers, sequences = parse_genbank(input_file)
        else:
            continue
        for geneid, translation in zip(headers, sequences):
            try:
                database.add_gene(geneid, genome, translation)
            except ValueError as err:
                print(f"Error {err.args[0]}: {err.args[1]}")

    return database


def dump_nr_translations(database, filepath):
    """
    Dump the non-redundant translations to a FASTA file for clustering.

    :param database: a database containing genes
    :type database: SequenceDB or Pham
    :param filepath: file to write non-redundant translations to
    :type filepath: pathlib.Path
    """
    headers, sequences = list(), list()
    for geneid, translation in database.nr_iterator:
        headers.append(geneid)
        sequences.append(translation)
    write_fasta(headers, sequences, filepath)


def merge_pre_and_hmm_phams(hmm_phams, pre_phams, database):
    """
    Merge first and second iteration phams.

    :param hmm_phams: clustered consensus sequences
    :type hmm_phams: dict
    :param pre_phams: clustered sequences (used to generate hmms)
    :type pre_phams: dict
    :param database: the database whose sequences were clustered
    :type database: SequenceDB
    :return: merged_phams
    """
    phams, lookup = list(), dict()
    for key, translations in pre_phams.items():
        pham = Pham()
        for translation in translations:
            lookup[translation] = key
            for geneid in database.get_translation_geneids(translation):
                genome = database.get_geneid_genome(geneid)
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
                   verbose, skip_hmm, tmp_dir, debug=False):
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
    :param debug: print verbose outputs to the console
    :type debug: bool
    :return: phams
    """
    if verbose:
        print("Creating MMseqs2 database...")

    nr_fasta = tmp_dir.joinpath("all_genes.fasta")
    iter_1_out = tmp_dir.joinpath("firstIter.txt")
    iter_2_out = tmp_dir.joinpath("secondIter.txt")

    dump_nr_translations(database, nr_fasta)
    mmseqs_createdb(nr_fasta, tmp_dir, debug=debug)

    if verbose:
        print("Performing sequence-sequence clustering...")

    clu_mode, sens, i, c, e = first_iter
    mmseqs_cluster(tmp_dir, clu_mode, sens, i, c, e, threads=threads,
                   debug=debug)

    if verbose:
        print("Parsing first iteration phams...")

    mmseqs_first_iter_cleanup(tmp_dir, iter_1_out, threads, debug=debug)

    first_iter_phams = parse_mmseqs_output(iter_1_out)

    if not skip_hmm:
        if verbose:
            print("Building HMMs from pre-phams...")

        mmseqs_result2profile(tmp_dir, threads, debug=debug)

        if verbose:
            print("Extracting consensus sequences from HMMs...")

        mmseqs_profile2consensus(tmp_dir, threads, debug=debug)

        if verbose:
            print("Performing consensus-HMM clustering...")

        clu_mode, sens, i, c, e = second_iter
        mmseqs_search(tmp_dir, i, c, e, threads=threads, debug=debug)

        mmseqs_clust(tmp_dir, clu_mode, threads, debug=debug)

        if verbose:
            print("Parsing second iteration phams...")

        mmseqs_second_iter_cleanup(tmp_dir, iter_2_out, threads, debug=debug)

        second_iter_phams = parse_mmseqs_output(iter_2_out)
    else:
        second_iter_phams = None

    return merge_pre_and_hmm_phams(second_iter_phams, first_iter_phams,
                                   database)


def write_pham_fastas(phams, fasta_dir, nr_only=False):
    """
    Dump pham sequences to FASTA files (one file per pham).

    :param phams: the assembled phams
    :type phams: list[Pham]
    :param fasta_dir: the directory where files should be written
    :type fasta_dir: pathlib.Path
    :param nr_only: only write non-redundant sequences to files
    :type nr_only: bool
    """
    for i, pham in enumerate(phams):
        fasta_path = fasta_dir.joinpath(f"pham_{i+1}.faa")
        if nr_only:
            dump_nr_translations(pham, fasta_path)
        else:
            headers, sequences = list(), list()
            for _, geneid, translation in pham:
                headers.append(geneid), sequences.append(translation)
            write_fasta(headers, sequences, fasta_path)


def map_pangenome(phams, n_genomes):
    """
    Map out which phams belong to the core, soft-core, shell, or cloud
    compartments of the pangenome.

    :param phams: the assembled gene phamilies
    :type phams: list[Pham]
    :param n_genomes: the number of genomes examined
    :type n_genomes: int
    :return: pangenome_map
    """
    pangenome_map = {"core": [], "soft-core": [], "shell": [], "cloud": []}

    for i, pham in enumerate(phams):
        phamid = f"pham_{i+1}"
        fraction = float(len(pham.genomes)) / n_genomes
        if fraction >= 0.99:
            pangenome_map["core"].append(phamid)
        elif fraction >= 0.95:
            pangenome_map["soft-core"].append(phamid)
        elif fraction >= 0.15:
            pangenome_map["shell"].append(phamid)
        else:
            pangenome_map["cloud"].append(phamid)

    return pangenome_map


def summarize(pangenome_map, phams, all_genomes, out_dir):
    """
    Create some of the Roary output files.

    :param pangenome_map: the core/soft/shell/cloud mappings
    :type pangenome_map: dict[str[list[Pham]]]
    :param phams: the assembled phams
    :type phams: list[Pham]
    :param all_genomes:
    :type all_genomes: list[str]
    :param out_dir: the directory where Roary-style files can go
    :type out_dir: pathlib.Path
    """
    summary_statistics = out_dir.joinpath("summary_statistics.txt")
    with open(summary_statistics, "w") as summary_writer:
        summary_writer.write(f"Core genes (99% <= strains <= 100%):\t"
                             f"{len(pangenome_map['core'])}\n")
        summary_writer.write(f"Soft core genes (95% <= strains < 99%):\t"
                             f"{len(pangenome_map['soft-core'])}\n")
        summary_writer.write(f"Shell genes (15% <= strains < 95%):\t"
                             f"{len(pangenome_map['shell'])}\n")
        summary_writer.write(f"Cloud genes (0% <= strains < 15%):\t"
                             f"{len(pangenome_map['cloud'])}\n")
        summary_writer.write(f"Total genes (0% <= strains <= 100%):\t"
                             f"{len(phams)}\n")

    core_genes = out_dir.joinpath("core_genes.txt")
    with open(core_genes, "w") as core_writer:
        for core_phamid in pangenome_map['core']:
            core_writer.write(f"{core_phamid}\n")

    soft_genes = out_dir.joinpath("soft_core_genes.txt")
    with open(soft_genes, "w") as soft_writer:
        for soft_phamid in pangenome_map['soft-core']:
            soft_writer.write(f"{soft_phamid}\n")

    shell_genes = out_dir.joinpath("shell_genes.txt")
    with open(shell_genes, "w") as shell_writer:
        for shell_phamid in pangenome_map['shell']:
            shell_writer.write(f"{shell_phamid}\n")

    cloud_genes = out_dir.joinpath("cloud_genes.txt")
    with open(cloud_genes, "w") as cloud_writer:
        for cloud_phamid in pangenome_map['cloud']:
            cloud_writer.write(f"{cloud_phamid}\n")

    strain_genes = out_dir.joinpath("strain_genes.tsv")
    with open(strain_genes, "w") as fh:
        for i, pham in enumerate(phams):
            phamid = f"pham_{i+1}"
            for genome in pham.genomes:
                for _ in pham.get_genome_geneids(genome):
                    fh.write(f"{phamid}\t{genome}\n")

    gene_presence_absence = out_dir.joinpath("gene_presence_absence.csv")
    with open(gene_presence_absence, "w") as roary_writer:
        header = HEADER
        header.extend(all_genomes)
        roary_writer.write(",".join(header)+"\n")
        for i, pham in enumerate(phams):
            phamid = f"pham_{i+1}"
            seq_per_iso = float(len(pham))/len(pham.genomes)
            output_row = [phamid, '', str(len(pham.genomes)),
                          str(len(pham)), str(seq_per_iso),
                          '', '', '', '', '',
                          str(pham.minimum_length()),
                          str(pham.maximum_length()),
                          str(pham.average_length())]
            for genome in all_genomes:
                if pham.in_genome(genome):
                    geneids = pham.get_genome_geneids(genome)
                    geneids = ";".join([x.split()[0] for x in geneids])
                    output_row.append(geneids)
                else:
                    output_row.append('')
            roary_writer.write(",".join(output_row) + "\n")
