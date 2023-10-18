def gc_content(seq: str):
    gc_total = 0
    for nucleotide in seq:
        if nucleotide.lower() == 'g' or nucleotide.lower() == 'c':
            gc_total += 1
    return gc_total*100/len(seq)


def is_in_gc_bounds(seq: str, gc_lower_bound: int, gc_upper_bound: int):
    if gc_lower_bound <= gc_content(seq) <= gc_upper_bound:
        return True
    else:
        return False


def is_in_length_bounds(seq: str, length_lower_bound: int, length_upper_bound: int):
    if length_lower_bound <= len(seq) <= length_upper_bound:
        return True
    else:
        return False


def is_above_quality_threshold(seq_quality, quality_threshold):
    total_q_score = 0
    for ascii_quality_sign in seq_quality:
        q_score = ord(ascii_quality_sign) - 33
        total_q_score += q_score
    if total_q_score/len(seq_quality) >= quality_threshold:
        return True
    else:
        return False
