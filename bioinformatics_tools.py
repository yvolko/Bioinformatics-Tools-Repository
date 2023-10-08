from typing import Dict, Tuple
from additional_modules.additional_modules_fastq_thresholding import gc_content, is_in_gc_bounds, is_in_length_bounds, is_above_quality_threshold
from additional_modules.additional_modules_protein_analysis import molecular_weight, one_letter_to_three, get_amino_acid_sum, codon_optimization, length, name_transform, is_amino_acid, brutto_count, is_length_divisible_by_3, is_amino_acid_three_letter, 
from additional_modules.additional_modules_run_dna_rna_tools import is_valid_dna_rna, transcribe, reverse, complement, reverse_complement, get_nucl_acid_type

def run_dna_rna_tools(*args: str):
    procedure = args[-1]
    seqs = args[:-1]

    for seq in seqs:
        is_valid_dna_rna(seq)

    processed_result = []
    for seq in seqs:
        if procedure == 'transcribe':
            processed_result.append(transcribe(seq))
        elif procedure == 'reverse':
            processed_result.append(reverse(seq))
        elif procedure == 'complement':
            processed_result.append(complement(seq))
        elif procedure == 'reverse_complement':
            processed_result.append(reverse_complement(seq))
        elif procedure == 'get_nucl_acid_type':
            processed_result.append(get_nucl_acid_type(seq))
        else:
            return 'Procedure is not defined in the function'

    if len(processed_result) == 1:
        return processed_result[0]
    return processed_result

def protein_analysis(
        *args: str, procedure: str, cell_type: str = None, letter_format: int = 1
) -> list:
    """
    Function protein_analysis:
    - calculates predicted molecular weight of amino acid sequences in kDa (procedure name: molecular_weight)
    - translate aa sequences from one-letter to three-letter code (procedure name: one_letter_to_three)
    - calculates total amount of each amino acid in the sequences (procedure name: get_amino_acid_sum)
    - makes DNA based codon optimization for the introduced amino acid sequences, support 3 types of cells:
      Esherichia coli, Pichia pastoris, Mouse (procedure name: codon_optimization)
    - calculates length of amino acid sequences (procedure name: length)
    - counts the number of atoms of each type in a sequence (procedure name: brutto_count)

    Arguments:
    - one or multiple string of protein sequences written one letter or three letter code (not mixed)
    - name of procedure as string
    - cell type (required only for codon_optimization procedure)
    - letter_format of code for the protein sequences as int: 1 for one letter, 3 for three letter code

    Return:
    - molecular_weight procedure returns list of floats
    - one_letter_to_three procedure returns list of strings
    - get_amino_acid_sum procedure returns list of dictionaries
    - codon_optimization procedure returns list of strings
    - length procedure returns list of int values
    - brutto_count procedure returns list of dictionaries with counts of atoms in the sequence
    """
    amino_acid_seqs = name_transform(args, letter_format)
    procedures = {
        "molecular_weight": molecular_weight,
        "one_letter_to_three": one_letter_to_three,
        "get_amino_acid_sum": get_amino_acid_sum,
        "codon_optimization": codon_optimization,
        "length": length,
        "brutto_count": brutto_count,
    }
    if procedure not in procedures.keys():
        raise ValueError("Requested procedure is not defined")
    elif procedure == "codon_optimization":
        return procedures.get(procedure)(amino_acid_seqs, cell_type)
    else:
        return procedures.get(procedure)(amino_acid_seqs)

    
def fastq_thresholding(seqs: Dict[str, Tuple[str, str]],
                       gc_bounds: Tuple[int, int] = (0, 100),
                       length_bounds: Tuple[int, int] = (0, 2**32),
                       quality_threshold: int = 0):
    thresholded_dict = {}
    if isinstance(gc_bounds, int):
        gc_lower_bound = 0
        gc_upper_bound = gc_bounds
    else:
        gc_lower_bound = gc_bounds[0]
        gc_upper_bound = gc_bounds[1]

    if isinstance(length_bounds, int):
        length_lower_bound = 0
        length_upper_bound = length_bounds
    else:
        length_lower_bound = length_bounds[0]
        length_upper_bound = length_bounds[1]

    for key, value in seqs.items():
        if is_in_gc_bounds(value[0], gc_lower_bound, gc_upper_bound) & \
                is_in_length_bounds(value[0], length_lower_bound, length_upper_bound) & \
                is_above_quality_threshold(value[1], quality_threshold):
            thresholded_dict[key] = value
    return thresholded_dict
