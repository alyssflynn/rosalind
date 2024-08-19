"""Collection of solutions for problems from rosalind.info.

https://rosalind.info/problems/list-view/
"""

import re
from collections.abc import Iterable

from rosalind.const import DNA_BASE_PAIR, FASTA_PAT, RNA_CODON_TABLE


def count_bases(dna: str) -> dict:
    """Returns the count of each base pair in the DNA or RNA string."""
    counter = {"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}
    for nuc in dna:
        counter[nuc] += 1

    return counter


def base_count_str(counts: dict) -> str:
    """Returns the base count as a string.

    NOTE: rosalind expects this format in one problem but it's not incredibly useful to include in this package.
    """
    return f"{counts['A']} {counts['C']} {counts['G']} {counts['T']}"


def transcribe(dna: str) -> str:
    """Returns the transcribed RNA string for a given DNA string (replaces T's with U's)."""
    return re.sub("T", "U", dna)


def complement(dna: str) -> str:
    """Returns the complementary base pair sequence for a DNA or RNA string."""
    return "".join(DNA_BASE_PAIR[nuc] for nuc in dna)


def reverse_complement(dna: str) -> str:
    """Returns the reverse complement of a DNA or RNA string."""
    return complement(dna[::-1])


def recurrence_relation(n: int, k: int) -> int:
    """Returns the total # of breeding pairs present after n months with each pair producing k offspring per generation.

    Begins with 1 non-breeding pair and in each generation, every pair of reproduction-age pairs produces a litter of k
    non-breeding pairs. Offspring take 1 month to mature to breeding age.

    https://rosalind.info/problems/fib/
    """
    adult_pairs = 0
    child_pairs = 1

    for _ in range(n):
        curr_adult_pairs = adult_pairs
        adult_pairs += child_pairs
        child_pairs = curr_adult_pairs * k

    return adult_pairs


def gc_content(dna: str) -> float:
    """Returns the GC content of a DNA string (% C or G)."""
    counts = count_bases(dna)
    return sum((counts["G"], counts["C"])) / sum(counts.values()) * 100


def clean_sequence(seq: str) -> str:
    """Removes non-DNA/RNA characters from string."""
    return re.sub(r"[^ACTGU]", "", seq)


def parse_fasta(fasta: str) -> dict[str, str]:
    """Parses string in FASTA format, returns dictionary of DNA sequences keyed on their label.

    https://rosalind.info/problems/gc/
    """
    return {label: clean_sequence(seq) for label, seq in re.findall(FASTA_PAT, fasta)}


def max_gc_fasta(fasta: str) -> tuple[str, float]:
    """Parses string in FASTA format, returns the label and GC content of the DNA string with the highest GC content.

    https://rosalind.info/problems/gc/
    """
    fasta_gc = {label: round(gc_content(seq), 5) for label, seq in parse_fasta(fasta).items()}
    label = max(fasta_gc.keys(), key=fasta_gc.get)
    return label, fasta_gc[label]


def count_point_mutations(s: str, t: str) -> int:
    """Returns the Hamming distance between s and t. Assumes s and t are of equal length.

    The Hamming distance is the number of corresponding symbols that differ in s and t.
    https://rosalind.info/problems/hamm/
    """
    hdist = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            hdist += 1
    return hdist


def probability_dominant(k: int, m: int, n: int) -> float:
    """Returns the probablility that an offspring with a dominant allele will spawn in a population.

    Given a population of k homozygous dominant, m heterozygous, and n homozygous recessive individuals, calculates the
    probability that two randomly selected mating organisms will produce an individual posessing a dominant allele.

    k: homozygous dominant individuals
    m: heterozygous individuals
    n: homozygous recessive individuals

    https://rosalind.info/problems/iprb/
    https://stackoverflow.com/questions/25119106/rosalind-mendels-first-law-iprb
    """
    # calculate the probability of recessive traits only
    total = k + m + n
    two_recess = (n / total) * ((n - 1) / (total - 1))
    two_hetero = (m / total) * ((m - 1) / (total - 1))
    het_recess = (n / total) * (m / (total - 1)) + (m / total) * (n / (total - 1))
    recess_prob = two_recess + two_hetero * 1 / 4 + het_recess * 1 / 2
    dom_prob = 1 - recess_prob  # take the complement
    return round(dom_prob, 5)


def translate_rna(rna: str) -> str:
    """Translates RNA sequence into amino acid protein sequence.

    https://rosalind.info/problems/prot/
    """
    rna = clean_sequence(rna)
    codons = ""
    for i in range(3, len(rna), 3):
        codons += RNA_CODON_TABLE[rna[i - 3 : i]]
    return codons


def find_motif(seq: str, motif: str) -> tuple[int]:
    """Returns the starting potions of each instance of `motif` in `seq`."""

    def scan(seq: str, motif: str) -> Iterable[int]:
        msize = len(motif)
        i = 0
        while i + msize <= len(seq):
            if seq[i : i + msize] == motif:
                yield i + 1
            i += 1

    return tuple(scan(seq, motif))


def consensus_profile(*seqs: str) -> tuple[str, dict]:
    """Returns a consensus string and profile matrix representing the frequency of each nucleotide at each position.

    A consensus string is a string formed by taking the most common symbol at each position.

    A profile matrix is a 4xn matrix P in which P1,j represents the number of times that 'A' occurs in the jth position
    of one of the strings, P2,j represents the number of times that C occurs in the jth position, and so on.
    """
    size = len(seqs[0])
    profile = {nuc: [0] * size for nuc in "ACGT"}
    consensus = ""

    for i in range(size):
        for seq in seqs:
            nuc = seq[i]
            profile[nuc][i] += 1

        counts = {nuc: profile[nuc][i] for nuc in profile}
        max_nuc = max(counts, key=counts.get)
        consensus += max_nuc

    return consensus, profile


def fasta_consensus_profile(fasta: str) -> tuple[str, dict]:
    """Given a collection of DNA strings in FASTA format, returns a consensus string and profile matrix."""
    data = parse_fasta(fasta)
    seqs = data.values()
    return consensus_profile(*seqs)
