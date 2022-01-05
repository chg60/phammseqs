"""A class derived from SequenceDB that is better suited for the kinds of
operations one might wish to do with smaller groups of related sequences."""

from phamerate.classes.database import SequenceDB


class Pham(SequenceDB):
    """
    A type of SequenceDB, with additional convenience methods for
    performing actions that are useful when considering relatively
    small numbers of related sequences, and not so much for whole
    meta- or pan-genomes.
    """
    def __init__(self):
        super().__init__()
        self._geneids_to_aligned_translations = dict()

    @property
    def alignment_iterator(self):
        """Alternative to __iter__(), which allows iteration over
        aligned sequences rather than the unaligned sequences.
        """
        alignments = self._geneids_to_aligned_translations.items()
        for geneid, translation in alignments:
            genome = self.get_geneid_genome(geneid)
            yield genome, geneid, translation

    def merge(self, other):
        """Merge phams from multiple rounds of clustering.

        :param other: another pham to be merged into this one
        :type other: Pham
        """
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        for genome, geneid, translation in other:
            try:
                self.add_gene(geneid, genome, translation)
            except ValueError:
                pass    # Tried adding a gene that was already in the pham

    def add_aligned_translation(self, aligned_translation):
        """Add an aligned translation for this database to keep track of.

        :param aligned_translation: an aligned translation
        :type aligned_translation: str
        :raise: KeyError
        """
        translation = aligned_translation.replace("-", "")
        geneids = self.get_translation_geneids(translation)
        for geneid in geneids:
            self._geneids_to_aligned_translations[geneid] = aligned_translation

    def in_genome(self, genome):
        """Check whether this phamily is present in the indicated genome.

        :param genome: the genome to be checked against
        :type genome: str
        return: in_genome
        """
        return genome in self.genomes

    def minimum_length(self, kind="nucleotide"):
        """Return the length of the shortest sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: min_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        min_length = min(lengths)

        if kind == "nucleotide":
            return min_length * 3 + 3

        return min_length

    def average_length(self, kind="nucleotide"):
        """Return the length of the average sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: avg_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        avg_length = float(sum(lengths))/len(lengths)

        if kind == "nucleotide":
            return avg_length * 3 + 3

        return avg_length

    def maximum_length(self, kind="nucleotide"):
        """Return the length of the longest sequence in the phamily. Can
        return either 'nucleotide' or 'translation' length.

        :param kind: the type of sequence length to return
        :type kind: str
        :return: max_length
        """
        if kind not in ("nucleotide", "translation"):
            raise ValueError(f"argument 'kind' must be one of ('nucleotide', "
                             f"'translation')")

        lengths = [len(x) for x in self.translations]
        max_length = max(lengths)

        if kind == "nucleotide":
            return max_length * 3 + 3

        return max_length

    def __lt__(self, other):
        if not isinstance(other, Pham):
            raise TypeError(f"cannot merge Pham with {type(other)}")

        return len(self) < len(other)
