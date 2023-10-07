def is_valid_dna_rna(seq: str) -> bool:
    set_of_seq = set(seq.upper())
    if set_of_seq > {'A', 'T', 'G', 'C', 'U'}:
        raise ValueError('Invalid alphabet')
    elif ('T' in set_of_seq) and ('U' in set_of_seq):
        raise ValueError('Invalid alphabet')
    else:
        return True


def transcribe(seq: str) -> str:
    transcribed_seq = ''
    for nucleotide in seq:
        if nucleotide == 'T':
            transcribed_seq += 'U'
        elif nucleotide == 't':
            transcribed_seq += 'u'
        else:
            transcribed_seq += nucleotide
    print(transcribed_seq)
    return transcribed_seq


def reverse(seq: str) -> str:
    reversed_seq = seq[::-1]
    return reversed_seq


def complement(seq: str) -> str:
    complement_seq = ''
    dna_complement_map = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g'
    }
    for nucleotide in seq:
        complement_seq += dna_complement_map[nucleotide]
    return complement_seq


def reverse_complement(seq: str) -> str:
    return complement(reverse(seq))


def get_nucl_acid_type(seq: str) -> str:
    set_of_seq = set(seq)
    if ('U' in set_of_seq) or ('u' in set_of_seq):
        return 'RNA'
    elif ('T' in set_of_seq) or ('t' in set_of_seq):
        return 'DNA'
    return 'ND'
