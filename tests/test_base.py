from rosalind import base_count_str, complement, count_bases, reverse_complement, transcribe


def test_count_bases():
    assert count_bases("GATTACA") == {"G": 1, "A": 3, "T": 2, "C": 1, "U": 0}


def test_base_count_str():
    assert base_count_str(count_bases("GATTACA")) == "3 1 1 2"


def test_transcribe():
    dna_count = count_bases("GATTACA")
    rna_count = count_bases(transcribe("GATTACA"))
    assert rna_count["U"] == dna_count["T"]
    assert rna_count["T"] == 0


def test_complement():
    assert complement("AAAACCCGGT") == "TTTTGGGCCA"


def test_reverse_complement():
    assert reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"
