import os
import re
import time
import requests
import sys
import pandas as pd
from typing import Tuple
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC
from additional_modules import transcribe, complement, molecular_weight
from io import StringIO
from dotenv import load_dotenv
from datetime import timedelta
from dataclasses import dataclass
from typing import List
from bs4 import BeautifulSoup


# DEFINE CLASSES OF BIOLOGICAL SEQUENCES
class BiologicalSequence(ABC):
    def __init__(self, sequence: str, alphabet: set):
        self.sequence = sequence
        self.alphabet = alphabet

    def __len__(self) -> int:
        """Returns the length of the sequence."""
        return len(self.sequence)

    def __getitem__(self, index) -> str:
        """Returns the item pointed by index or a slice of the sequence."""
        return self.sequence[index]

    def __str__(self) -> str:
        """Returns a string representation of the sequence."""
        return str(self.sequence)

    def check_alphabet(self) -> bool:
        """Checks if the sequence consists of allowed symbols."""
        return set(self.sequence) <= self.alphabet


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence, alphabet):
        super().__init__(sequence, alphabet)

    def gc_content(self) -> float:
        return gc_fraction(self.sequence)


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence, alphabet={'A', 'T', 'C', 'G'})

    def transcribe(self):
        return RNASequence(transcribe(self.sequence))

    def complement(self):
        return DNASequence(complement(self.sequence))


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence, alphabet={'A', 'C', 'G', 'U'})


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        super().__init__(sequence, alphabet=set("ACDEFGHIKLMNPQRSTVWY"))

    def molecular_weight(self):
        return molecular_weight(self.sequence)


# DEFINE CLASSES FOR GENSCAN
@dataclass
class GenomicElement:
    """
        Represents a genomic element (intron or exon) with
        number, start, and end position.

        Attributes:
            number (int): The identifier number of the genomic element,
            represents order of exons or introns.
            start (int): The starting position of the genomic element.
            end (int): The ending position of the genomic element.
    """
    number: int
    start: int
    end: int


@dataclass
class GenscanOutput:
    """
        Represents the output of Genscan prediction, including status,
        lists of coding sequences (CDS), introns, and exons.

        Attributes:
            status (str): Information about query to Genscan.
            cds_list (List[str]): A list of coding sequences (CDS)
            predicted by Genscan.
            intron_list (List[GenomicElement]): A list of introns.
            exon_list (List[GenomicElement]): A list of exons.
    """
    status: str
    cds_list: List[str]
    intron_list: List[GenomicElement]
    exon_list: List[GenomicElement]


# FASTA FILTRATOR
def filter_fastq(input_path: str, output_filename: str = '',
                 gc_bounds: Tuple[int, int] = (0, 100),
                 length_bounds: Tuple[int, int] = (0, 2 ** 32),
                 quality_threshold: int = 0):
    """
        filter_fastq(input_path: str, output_filename: str = '',
                 gc_bounds: Tuple[int, int] = (0, 100),
                 length_bounds: Tuple[int, int] = (0, 2 ** 32),
                 quality_threshold: int = 0)
        :param input_path: path to fastq file
        :param output_filename: desired name for the filtered fasta file
        (if not given _output is added to input file name)
        :param gc_bounds: borders of GC-content that will be used
        to filter sequence
        :param length_bounds: borders of sequence length that will
        be used to filter sequence
        :param quality_threshold: borders of quality that will
        be used to filter sequence (mean quality of the sequence
        is considered to pass the threshold)
        All the borders include upper and lower values.

        :return: dictionary of the same structure as input,
        that only contains sequences that passed all the
        thresholding
        """

    # Prepare output directory and output file name if not given

    dir_name = os.path.dirname(input_path)

    if output_filename == '':
        output_filename = re.sub(r'\.[^.]*$', '',
                                 os.path.basename(input_path)) + "_output"

    if not os.path.exists(os.path.join(dir_name, 'fastq_filtrator_resuls')):
        os.makedirs(os.path.join(dir_name, 'fastq_filtrator_resuls'))

    # Make sure all the filers are in the correct format

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

    # Read fasta recored, apply filters and write to the output file

    with open(input_path, "r") as fasta_file:
        with open(os.path.join(dir_name, 'fastq_filtrator_resuls',
                               f'{output_filename}.fasta'), mode='x') as file:
            for record in SeqIO.parse(fasta_file, "fastq"):
                if is_in_gc_bounds(record, gc_lower_bound, gc_upper_bound) \
                    and is_in_length_bounds(record, length_lower_bound,
                                            length_upper_bound) and \
                        is_above_quality_threshold(record, quality_threshold):
                    SeqIO.write(record, file, "fasta")


def is_in_gc_bounds(seq, gc_lower_bound: int, gc_upper_bound: int) -> bool:
    if gc_lower_bound <= gc_fraction(seq.seq) * 100 <= gc_upper_bound:
        return True
    return False


def is_in_length_bounds(seq, length_lower_bound: int, length_upper_bound: int):
    if length_lower_bound <= len(seq) <= length_upper_bound:
        return True
    return False


def is_above_quality_threshold(seq, quality_threshold):
    total_q_score = sum(seq.letter_annotations["phred_quality"])
    if total_q_score/len(seq) >= quality_threshold:
        return True
    return False


# TELEGRAM LOGGER (decorator for the functions to send logs with the chatbot)
dotenv = load_dotenv("keys.env")
token = os.getenv("TG_API_TOKEN")


def send_message(msg: str, chat_id: str, token: str) -> None:
    """Sends message to the telegram chatbot with given chat id"""
    send_msg = f"https://api.telegram.org/{token}/sendMessage?" \
               f"chat_id={chat_id}&text={msg}&parse_mode=HTML"
    requests.get(send_msg)


def send_file(func_name: str, file, chat_id: str, token: str) -> None:
    """Sends file to the telegram chatbot with given chat id"""
    file_content = file.getvalue()
    files = {'document': (f"{func_name}.log", file_content)}
    file_to_send = f"https://api.telegram.org/{token}/sendDocument?" \
                   f"chat_id={chat_id}&document="
    requests.post(file_to_send, files=files, timeout=2.50)


def telegram_logger(chat_id: str):
    """
    Decorator for the function. Sends to the telegram chatbot log file of
    the function execution and measures and sends time of the function run.
    :param chat_id: id of chatbot, string
    :return: wrapped function
    """

    def wrapper(func):
        def wrapped(*args, **kwargs):
            had_exception = False
            stdout_backup = sys.stdout
            output_buffer = StringIO()
            sys.stdout = output_buffer
            stderr_backup = sys.stderr
            sys.stderr = output_buffer
            start_time = time.time()
            try:
                result = func(*args, **kwargs)
            except Exception as e:
                had_exception = True
                exception_name = type(e).__name__
                exception_text = str(e)
            end_time = time.time()
            output = output_buffer.getvalue()
            sys.stdout = stdout_backup
            sys.stderr = stderr_backup
            file_buffer = StringIO()
            file_buffer.write(output)
            file_buffer.seek(0)

            if file_buffer.getvalue():
                send_file(func.__name__, file_buffer, chat_id, token)
            if had_exception:
                send_message(f"😔Function <code>{func.__name__}</code> failed "
                             f"with an exception: <code>{exception_name}: "
                             f"{exception_text}</code>", chat_id,
                             token)
                exit()
            run_time = end_time - start_time
            if run_time < 86400:
                run_time = str(timedelta(seconds=run_time))
            else:
                run_time = str(timedelta(seconds=int(run_time)))
            send_message(f"🎉Function <code>{func.__name__}</code> "
                         f"successfully finished "
                         f"in <code>{run_time}</code>", chat_id, token)
            return result

        return wrapped

    return wrapper


# GENSCAN - identify complete gene structures in genomic DNA
# (predicts the location of genes and their exon-intron boundaries
# in genomic sequences from a variety of organisms)
def run_genscan(sequence: str = None, sequence_file=None,
                organism: str = "Vertebrate", exon_cutoff: float = 1.00,
                sequence_name: str = ""):
    """
    Sends post query to Genescan and parse the response.

    :param sequence: nucleotide sequence
    :param sequence_file: file with nucleotide sequence
    :param organism: one of "Vertebrate", "Arabidopsis", "Maize"
    :param exon_cutoff: suboptimal exon cutoff, one of (1.0, 0.5, 0.25, 0.1,
    0.05, 0.02, 0.01)
    :param sequence_name: name of the sequence
    :return: GenscanOutput objet
    """

    if sequence_file is not None:
        with open(sequence_file, 'r') as file:
            sequence = file.read()
        sequence = [line.rstrip('\n') for line in sequence]
        sequence = ''.join(sequence)

    files = {
        '-o': (None, organism),
        '-e': (None, exon_cutoff),
        '-n': (None, sequence_name),
        '-p': (None, 'Predicted peptides only'),
        '-s': (None, sequence),
    }
    response = requests.post('http://hollywood.mit.edu/cgi-bin/'
                             'genscanw_py.cgi', data=files)

    soup = BeautifulSoup(response.content, "html.parser")
    predictions = soup.find_all("pre")
    predictions = str(predictions).splitlines()
    predictions = [line for line in predictions if line.strip()]

    # Get status
    status = predictions[1]
    predictions = predictions[7:]

    # Get exon table
    exons = []
    for line in predictions:
        if line.startswith('Suboptimal'):
            break
        exons.append(line)
    split_exons = [line.split() for line in exons if ('Init' in line or
                                                      'Intr' in line or
                                                      'Term' in line)]

    exon_df = pd.DataFrame(split_exons, columns=["Gn.Ex", "Type", "S",
                                                 "Begin", "End", "Len",
                                                 "Fr", "Ph", "I/Ac", "Do/T",
                                                 "CodRg", "P", "Tscr"])

    def get_exons(row):
        num = int(row['Gn.Ex'].split(".")[1])
        start = int(row['Begin'])
        end = int(row['End'])

        return GenomicElement(num, start, end)

    all_exons = list(exon_df.apply(get_exons, axis=1))

    introns = list()
    for i in range(1, len(all_exons)):
        end = all_exons[i].start
        start = all_exons[i - 1].end
        introns.append(GenomicElement(i, start, end))

    # Get protein predictions
    cut_ind = 0
    for ind, line in enumerate(predictions):
        if line.startswith('Predicted peptide sequence(s):'):
            cut_ind = ind + 2
    predictions = predictions[cut_ind:-1]
    protein_prediction = ''.join(predictions)
    if '&' in protein_prediction:
        protein_prediction = protein_prediction.split('&')
        for ind, pp in enumerate(protein_prediction):
            if '_aa' in pp:
                protein_prediction[ind] = pp.split('_aa')[1]

    return GenscanOutput(status, protein_prediction, introns, all_exons)


sequence = 'DYEVTFTEDKINAL'
peptide = AminoAcidSequence(sequence)
print(peptide.molecular_weight())