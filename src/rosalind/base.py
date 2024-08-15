"""Blah.

Yadda yadda.
"""

import re

DNA_BASE_PAIR = {"A": "T", "C": "G", "G": "C", "T": "A"}
RNA_BASE_PAIR = {"U": "T", "C": "G", "G": "C", "T": "U"}


def count_bases(dna: str) -> dict:
    """TODO."""
    counter = {"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}
    for nuc in dna:
        counter[nuc] += 1

    return counter


def base_count_str(counts: dict) -> str:
    """TODO."""
    return f"{counts['A']} {counts['C']} {counts['G']} {counts['T']}"


def transcribe(dna: str) -> str:
    """TODO."""
    return re.sub("T", "U", dna)


def complement(dna: str) -> str:
    """TODO."""
    return "".join(DNA_BASE_PAIR[nuc] for nuc in dna)


def reverse_complement(dna: str) -> str:
    """TODO."""
    return complement(dna[::-1])
