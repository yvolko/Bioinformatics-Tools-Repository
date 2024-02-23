amino_short_names_dic = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "V": "Val",
    "H": "His",
    "G": "Gly",
    "Q": "Gln",
    "E": "Glu",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "P": "Pro",
    "S": "Ser",
    "Y": "Tyr",
    "T": "Thr",
    "W": "Trp",
    "F": "Phe",
    "C": "Cys",
}

amino_names_dic = {
    "ala": "A",
    "arg": "R",
    "asn": "N",
    "asp": "D",
    "val": "V",
    "his": "H",
    "gly": "G",
    "gln": "Q",
    "glu": "E",
    "ile": "I",
    "leu": "L",
    "lys": "K",
    "met": "M",
    "pro": "P",
    "ser": "S",
    "tyr": "Y",
    "thr": "T",
    "trp": "W",
    "phe": "F",
    "cys": "C",
}

amino_names_dic_reverse = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Val": "V",
    "His": "H",
    "Gly": "G",
    "Gln": "Q",
    "Glu": "E",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Pro": "P",
    "Ser": "S",
    "Tyr": "Y",
    "Thr": "T",
    "Trp": "W",
    "Phe": "F",
    "Cys": "C",
}

amino_weights = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.16,
    "E": 147.13,
    "Q": 146.15,
    "G": 75.07,
    "H": 155.16,
    "I": 131.18,
    "L": 131.18,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

amino_brutto = {
    "A": (3, 7, 1, 2, 0),
    "R": (6, 14, 4, 2, 0),
    "N": (4, 8, 2, 3, 0),
    "D": (4, 7, 1, 4, 0),
    "V": (5, 11, 1, 2, 0),
    "H": (6, 9, 3, 2, 0),
    "G": (2, 5, 1, 2, 0),
    "Q": (5, 10, 2, 3, 0),
    "E": (5, 9, 1, 4, 0),
    "I": (6, 13, 1, 2, 0),
    "L": (6, 13, 1, 2, 0),
    "K": (6, 14, 2, 2, 0),
    "M": (5, 11, 1, 2, 1),
    "P": (5, 9, 1, 2, 0),
    "S": (3, 7, 1, 3, 0),
    "Y": (9, 11, 1, 3, 0),
    "T": (4, 9, 11, 1, 3, 0),
    "W": (11, 12, 2, 2, 0),
    "F": (9, 11, 1, 2, 0),
    "C": (3, 7, 1, 2, 1),
}

ecoli_triplets = {
    "A": "GCG",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGC",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCG",
    "Q": "CAG",
    "R": "CGT",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAT",
}

ppastoris_triplets = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAG",
    "L": "TTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCA",
    "Q": "CAA",
    "R": "AGA",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAC",
}

mouse_triplets = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "CGG",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
}


def molecular_weight(amino_acid_seqs: list) -> list:
    """
    Calculates predicated molecular weight of aa sequences.

     Arguments:
    - amino_acid_seqs (list): list of string with the protein sequences

    Return:
    - List of floats corresponding to the molecular weight in kDa
    """
    molecular_weights = []
    for seq in amino_acid_seqs:
        total_weight = 0
        for aa in seq:
            aa = aa.upper()
            total_weight += amino_weights[aa]
        molecular_weights.append(round(total_weight / 1000, 2))
    return molecular_weights


def one_letter_to_three(amino_acid_seqs: list) -> list:
    """
    Translates one-letter coded amino acid sequences to three-letter coded
    Arguments:
    - amino_acid_seqs (list): list of string with the protein sequences

    Return:
    - List of of strings with three-letter coded sequences
    """
    three_letters_seqs = []
    for seq in amino_acid_seqs:
        three_letters_seq = []
        for aa in seq:
            aa = aa.upper()
            three_letters_seq.append(amino_short_names_dic[aa])
        three_letters_seqs.append("".join(three_letters_seq))
    return three_letters_seqs


def get_amino_acid_sum(protein_sequences: list) -> list:
    """
    Counts the amount of each amino acid in the injected protein sequences

    Arguments:
    - protein_sequences (list): list of injected protein sequence

    Return:
    - List of dictionary with amino acid amount"""
    result = []
    for protein_sequence in range(len(protein_sequences)):
        amino_acid_count = {
            "A": 0,
            "C": 0,
            "D": 0,
            "E": 0,
            "F": 0,
            "G": 0,
            "H": 0,
            "I": 0,
            "K": 0,
            "L": 0,
            "M": 0,
            "N": 0,
            "P": 0,
            "Q": 0,
            "R": 0,
            "S": 0,
            "T": 0,
            "V": 0,
            "W": 0,
            "Y": 0,
        }
        for amino_acid in protein_sequences[protein_sequence]:
            amino_acid_count[amino_acid] += 1
        result.append(amino_acid_count)
    return result


def codon_optimization(protein_sequences: list, cell_type: str) -> list:
    """
    Makes codon-optimized DNA based on the introduced amino acid sequences for 3 types of cells:
    Esherichia coli, Pichia pastoris, Mouse

    Arguments:
    - protein_sequences (list): list of injected protein sequence
    - cell_type (str): user-entered cell type for codon optimization

    Return:
    - List of codon-optimized DNA"""

    if cell_type == "Esherichia coli" or cell_type == "E.coli":
        codon_optimization_ecoli = []
        replacer_ecoli = ecoli_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_ecoli += [
                "".join([replacer_ecoli(n, n) for n in protein_sequences[amino_acid]])
            ]
        return codon_optimization_ecoli

    if cell_type == "Pichia pastoris" or cell_type == "P.pastoris":
        codon_optimization_ppastoris = []
        replacer_ppastoris = ppastoris_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_ppastoris += [
                "".join(
                    [replacer_ppastoris(n, n) for n in protein_sequences[amino_acid]]
                )
            ]
        return codon_optimization_ppastoris

    if cell_type == "Mouse" or cell_type == "mouse":
        codon_optimization_mouse = []
        replacer_mouse = mouse_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_mouse += [
                "".join([replacer_mouse(n, n) for n in protein_sequences[amino_acid]])
            ]
        return codon_optimization_mouse
    else:
        raise ValueError(
            f'Type {cell_type} is not supported. The following types of organisms are available for codon optimization: Esherichia coli, Pichia pastoris, Mouse'
        )


def length(seqs: list) -> list:
    """
    Counts total length of amino acid sequence.

    Arguments:
    - seqs (list): list of string with the protein sequences

    Return:
    - list of int values corresponding to the length of sequences"""
    result = [len(seq) for seq in seqs]
    return result


def name_transform(seqs: tuple, letter_format: int) -> list:
    """
    Transforms the amino acid sequences given to protein_analysis function from three-letter code to one-letter code,
    makes sequences unified (for one-letter letter_format all letters to upper and
    for three-letter letter_format to lower).

    Arguments:
      - seqs (tuple): tuple of string with the protein sequences

      Return:
      - list of strings with the transformed sequences"""
    result = []
    multiple_of_three = []
    test_three_letters = []
    if letter_format == 1:
        for seq in seqs:
            multiple_of_three.append(is_length_divisible_by_3(seq))
            test_three_letters.append(is_amino_acid_three_letter(seq))
            seq = seq.upper()
            for letter in seq:
                if is_amino_acid(letter):
                    pass
            result.append(seq)
        if all(multiple_of_three) and all(test_three_letters):
            print(
                "Warning: all your sequences are similar to three-letter ones. Check the letter_format value"
            )
        return result
    elif letter_format == 3:
        for seq in seqs:
            seq = seq.lower()
            seq3 = [seq[i: i + 3] for i in range(0, len(seq), 3)]
            for triplet in seq3:
                if is_amino_acid(triplet):
                    pass
            seq_transformed = "".join([amino_names_dic.get(seq) for seq in seq3])
            result.append(seq_transformed)
        return result
    else:
        raise ValueError(
            "Error unsupported letter_format. Only letter_formats 1 and 3 are supported"
        )


def is_amino_acid(input_amino: str) -> bool:
    """
    Checks whether the entered string is an amino acid (either three-letter encoding or one-letter encoded).

    Arguments:
      - input_amino (str): string corresponding to one amino acid (in three-letter code or one-letter code)

      Return:
      - bool: True if amino acid is a valid amino acid, otherwise ValueError is amino acid is not correct
    """
    if len(input_amino) == 1:
        letter = input_amino
        if letter not in amino_short_names_dic.keys():
            raise ValueError(f"Error {letter} is not an amino acid. Correct your input")
        return True
    elif len(input_amino) == 3:
        triplet = input_amino
        if triplet not in amino_names_dic.keys():
            raise ValueError(
                f"Error {triplet} is not an amino acid. Correct your input"
            )
        return True
    else:
        raise ValueError(
            f"Error {input_amino} is incorrect form of amino acid notation. Correct your input"
        )


def brutto_count(seqs: list) -> list:
    """
    Calculates the brutto formula of the amino acid sequences.

    Arguments:
      - seqs (list): list of string with the protein sequences

      Return:
      - list of dictionaries with counts of each elemet included (elements C,H,N,O,S)"""
    elements = ["C", "H", "N", "O", "S"]
    result = []
    for seq in seqs:
        brutto_list = [amino_brutto.get(letter) for letter in seq]
        brutto_pair = list(zip(*brutto_list))
        brutto = [sum(i) for i in brutto_pair]
        brutto_dict = dict(zip(elements, brutto))
        result.append(brutto_dict)
    return result


def is_length_divisible_by_3(seq: str) -> bool:
    """
    Checks if the sequence is divisible by three.

    Arguments:
      - seq (str): string of protein sequence

      Return:
      - bool: True if sequence is divisible by three, otherwise False"""
    seq_len = len(seq)
    if seq_len % 3 == 0:
        return True
    else:
        return False


def is_amino_acid_three_letter(seq: str) -> bool:
    """
    Checks whether all elements of a sequence are three-letter amino acid symbols.

    Arguments:
      - seq (str): string of protein sequence

    Return:
      - bool: True if sequence is corresponding to the valid three-letter amino acid, otherwise False
    """
    seq = seq.lower()
    seq3 = [seq[i: i + 3] for i in range(0, len(seq), 3)]
    for triplet in seq3:
        if triplet not in amino_names_dic.keys():
            return False
        else:
            return True
